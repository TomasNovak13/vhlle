class EoS;
class EoSaux;
class EoSChiral : public EoS {
private:
 EoSaux *eossmall, *eosbig;

public:
  double psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8, psi9, psi10;
 EoSChiral(void);
 ~EoSChiral(void);

 virtual void eos(double e, double nb, double nq, double ns, double &_T,
                  double &_mub, double &_muq, double &_mus, double &_p);
 virtual double p(double e, double nb, double nq, double ns);
};
