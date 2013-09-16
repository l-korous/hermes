/// Class for Heating induced vortex, linear initial condition.
class InitialSolutionShockTube : public ExactSolutionScalar<double>
{
public:
  InitialSolutionShockTube(MeshSharedPtr mesh, int component) : ExactSolutionScalar<double>(mesh), component(component)
  {
  };

  virtual double value (double x, double y) const 
  {
    double x_c;
    this->element->get_center(x_c, y);
    if(std::abs(x - 0.5) < 1e-6)
    {
      if(x_c < 0.4999)
        x = 0.49;
      else
        x = 0.51;
    }

    if(x < 0.5)
    {
      switch(this->component)
      {
      case 0:
        return 1.0;
      case 1:
      case 2:
        return 0.;
      case 3:
        return QuantityCalculator::calc_energy(1.0, 0., 0., 1.0, KAPPA);
      };
    }
    else
    {
      switch(this->component)
      {
      case 0:
        return 0.125;
      case 1:
      case 2:
        return 0.;
      case 3:
        return QuantityCalculator::calc_energy(0.125, 0., 0., 0.1, KAPPA);
      };
    }
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const
  {
    dx = 0.;
    dy = 0.;
  };

  virtual Ord ord(double x, double y) const
  {
    return Ord(1);
  }

  MeshFunction<double>* clone() const
  {
    return new InitialSolutionShockTube(this->mesh, this->component);
  }

  // Value.
  int component;
};