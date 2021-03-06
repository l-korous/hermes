// This file is part of Hermes2D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "api.h"
#include "config.h"

#include "callstack.h"
#include "common.h"
#include "exceptions.h"
#include "matrix.h"
#include "solvers/interfaces/paralution_solver.h"
#if defined __GNUC__ && defined HAVE_BFD
#include <signal.h>
#include "third_party/backtrace.c"
#endif
namespace Hermes
{
#if defined __GNUC__ && defined HAVE_BFD
  void hdl(int sig, siginfo_t *siginfo, void *context)
  {
    printf("Signal handler exception caught, Signal [number, error_number, code]: [%i, %i, %i].\n.", siginfo->si_signo, siginfo->si_errno, siginfo->si_code);

    void *array[100];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 100);

    // print out all the frames to stderr
    Hermes::backtrace_symbols_fd(array, size, 0);

    throw Hermes::Exceptions::Exception("signal caught");
  }
#endif

  Api::Parameter::Parameter(int default_val)
  {
    this->default_val = default_val;
    this->user_set = false;
  }

  Api::Api()
  {
#if defined __GNUC__ && defined HAVE_BFD
    act.sa_sigaction = &hdl;
    act.sa_flags = SA_SIGINFO;
    sigaction(SIGSEGV, &this->act, 0);
#endif

    // Insert parameters.
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*>(Hermes::numThreads, new Parameter(NUM_THREADS)));
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*>(Hermes::matrixSolverType, new Parameter(SOLVER_UMFPACK)));
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*>(Hermes::directMatrixSolverType, new Parameter(SOLVER_UMFPACK)));
#ifdef _DEBUG
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*>(Hermes::showInternalWarnings, new Parameter(1)));
#else
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*> (Hermes::showInternalWarnings, new Parameter(0)));
#endif
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*>(Hermes::useAccelerators, new Parameter(1)));
    this->parameters.insert(std::pair<HermesCommonApiParam, Parameter*>(Hermes::checkMeshesOnLoad, new Parameter(1)));

    // Set handlers.
#ifdef WITH_PARALUTION
    this->setter_handlers.insert(std::pair<HermesCommonApiParam, typename Api::SetterHandler>(Hermes::numThreads, &ParalutionInitialization::set_threads_paralution));
    this->value_setter_handlers.insert(std::pair<std::pair<HermesCommonApiParam, int>, typename Api::SetterHandler>(std::pair<HermesCommonApiParam, int>(Hermes::matrixSolverType, SOLVER_PARALUTION_ITERATIVE), &ParalutionInitialization::init_paralution));
    this->value_setter_handlers.insert(std::pair<std::pair<HermesCommonApiParam, int>, typename Api::SetterHandler>(std::pair<HermesCommonApiParam, int>(Hermes::matrixSolverType, SOLVER_PARALUTION_AMG), &ParalutionInitialization::init_paralution));
    this->change_handlers.insert(std::pair<std::pair<HermesCommonApiParam, int>, typename Api::SetterHandler>(std::pair<HermesCommonApiParam, int>(Hermes::matrixSolverType, SOLVER_PARALUTION_ITERATIVE), &ParalutionInitialization::deinit_paralution));
    this->change_handlers.insert(std::pair<std::pair<HermesCommonApiParam, int>, typename Api::SetterHandler>(std::pair<HermesCommonApiParam, int>(Hermes::matrixSolverType, SOLVER_PARALUTION_AMG), &ParalutionInitialization::deinit_paralution));
#endif

    // Initialize TCMalloc (this line also serves for TCMalloc not to be linker-optimized out).
#ifdef WITH_TC_MALLOC
    ::tc_set_new_mode(1);
#endif
  }

  Api::~Api()
  {
    // Call the change_handlers upon exit.
    for (std::map<HermesCommonApiParam, Parameter*>::iterator param_iterator = this->parameters.begin(); param_iterator != this->parameters.end(); param_iterator++)
    {
      std::map<std::pair<HermesCommonApiParam, int>, SetterHandler>::iterator change_handler_iterator = this->change_handlers.find(std::pair<HermesCommonApiParam, int>(param_iterator->first, param_iterator->second->user_set ? param_iterator->second->user_val : param_iterator->second->default_val));
      if (change_handler_iterator != this->change_handlers.end())
        change_handler_iterator->second();
    }

    for (std::map<HermesCommonApiParam, Parameter*>::const_iterator it = this->parameters.begin(); it != this->parameters.end(); ++it)
      delete it->second;
  }

  int Api::get_integral_param_value(HermesCommonApiParam param)
  {
    if (this->parameters.find(param) == parameters.end())
      throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);
    if (this->parameters.find(param)->second->user_set)
      return this->parameters.find(param)->second->user_val;
    else
      return this->parameters.find(param)->second->default_val;
  }

  void Api::set_integral_param_value(HermesCommonApiParam param, int value)
  {
    if (this->parameters.find(param) == parameters.end())
      throw Hermes::Exceptions::Exception("Wrong Hermes::Api parameter name:%i", param);

    // If already set, we might have to call some custom handler.
    if (this->parameters.find(param)->second->user_set)
    {
      std::map<std::pair<HermesCommonApiParam, int>, SetterHandler>::iterator change_handler_iterator = this->change_handlers.find(std::pair<HermesCommonApiParam, int>(param, this->parameters.find(param)->second->user_val));
      if (change_handler_iterator != this->change_handlers.end())
      {
        change_handler_iterator->second();
      }
    }
    this->parameters.find(param)->second->user_set = true;
    this->parameters.find(param)->second->user_val = value;

    // And now we might have to call some handler of the new_ value setting.
    // First - generic.
    std::map<HermesCommonApiParam, SetterHandler>::iterator setter_handler_iterator = this->setter_handlers.find(param);
    if (setter_handler_iterator != this->setter_handlers.end())
    {
      setter_handler_iterator->second();
    }

    // Second - specific.
    std::map<std::pair<HermesCommonApiParam, int>, SetterHandler>::iterator value_setter_handler_iterator = this->value_setter_handlers.find(std::pair<HermesCommonApiParam, int>(param, this->parameters.find(param)->second->user_val));
    if (value_setter_handler_iterator != this->value_setter_handlers.end())
    {
      value_setter_handler_iterator->second();
    }
  }

#if defined(WIN32) || defined(_WINDOWS)
  __declspec(dllexport) Hermes::Api HermesCommonApi;
#else
  Hermes::Api HermesCommonApi;
#endif
}
