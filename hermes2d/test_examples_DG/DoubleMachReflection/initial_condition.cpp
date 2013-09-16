class CustomInitialCondition : public ExactSolutionScalar<double>
{
public:
  CustomInitialCondition(MeshSharedPtr mesh, int component, double kappa) : ExactSolutionScalar<double>(mesh), component(component), kappa(kappa)
  {
this->time = 0.;
    sqrt3 = std::sqrt(3.0);
    first_vel = 8.25*std::cos(M_PI/6.)*8.0;
    second_vel = -8.25*std::sin(M_PI/6.)*8.0;
    first_e = QuantityCalculator::calc_energy(8.0, first_vel,second_vel, 116.5, kappa);
    second_e = QuantityCalculator::calc_energy(1.4, 0.0 ,0.0, 1.0, kappa);
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const 
  {
    dx = 0.0;
    dy = 0.0;
  }

  virtual double value (double x, double y) const
  {
    switch(this->component)
    {
    case 0:
      if(x < (1./6.) + (y * (1 + 20 * this->time) / sqrt3))
        return 8.0;
      else		
        return 1.4;
      break;
    case 1:
      if(x < (1./6.) + (y * (1 + 20 * this->time) / sqrt3))
        return first_vel;
      else		
        return 0.0;
      break;
    case 2:
      if(x < (1./6.) + (y * (1 + 20 * this->time) / sqrt3))
        return second_vel;
      else		
        return 0.0;
      break;
    case 3:
      if(x < (1./6.) + (y * (1 + 20 * this->time) / sqrt3))
        return first_e;
      else		
        return second_e;
      break;
    }
  }

  virtual Ord ord(double x, double y) const
  {
    return Ord(2);
  }

  MeshFunction<double>* clone() const
  {
    return new CustomInitialCondition(this->mesh, this->component, this->kappa);
  }

  int component;
  double time;
  double sqrt3, first_vel, second_vel, first_e, second_e;
  double kappa;
};
