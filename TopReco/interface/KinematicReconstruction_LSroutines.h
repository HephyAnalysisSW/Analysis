#ifndef KinematicReconstruction_LSroutines_h
#define KinematicReconstruction_LSroutines_h

#include <vector>

#include <TLorentzVector.h>

class TH1F;
class TF1;
class TF2;

#include "classesFwd.h"





class KinematicReconstruction_LSroutines{

public:
    KinematicReconstruction_LSroutines();
    KinematicReconstruction_LSroutines(const float& mass_Wp, const float& mass_Wm);
    KinematicReconstruction_LSroutines(const float& mass_top, const float& mass_topbar, 
                                       const float& mass_b, const float& mass_bbar, 
                                       const float& mass_w, const float& mass_wbar, 
                                       const float& mass_al, const float& mass_l);
    ~KinematicReconstruction_LSroutines();
    
    void fDelete()const;
    
    
    void ini(const float& mass_Wp, const float& mass_Wm);
    
    void setConstraints(const TLorentzVector& LV_al, 
                        const TLorentzVector& LV_l, 
                        const TLorentzVector& LV_b, 
                        const TLorentzVector& LV_bbar,
                        const float& missPx,
                        const float& missPy
                       );
    
    void setConstraints(const LV& LV_al, 
                        const LV& LV_l, 
                        const LV& LV_b, 
                        const LV& LV_bbar,
                        const float& missPx,
                        const float& missPy
                       );
    
    int getNsol()const;
    
    struct TopSolution{
        float dR;
        float dN;
        TLorentzVector top;
        TLorentzVector topbar;
        TLorentzVector neutrino;
        TLorentzVector neutrinobar;
        TLorentzVector wp;
        TLorentzVector wm;
        float x1;
        float x2;
        float mtt;
        float weight;
        
    };
    const std::vector<TopSolution>* getTtSol()const;
    
    void setTrueInfo(const LV& LV_Top,const LV& LV_AntiTop,const LV& LV_Neutrino, const LV& LV_AntiNeutrino);
    void sortBy(std::string ch);
    void print()const;
    
private:
    void filldR();
    void filldN();
    void swapTopSol(TopSolution& sol1, TopSolution& sol2)const;
    void sortTopSol(std::vector<TopSolution>& v)const;
    void doAll();
    void topRec(const float& px_neutrino);
    void findCoeff(float* const koeficienty);
    void quartic_equation(const float& h0, const float& h1, const float& h2, const float& h3, const float& h4, std::vector<float>& v)const;
    void cubic_equation(const float& a, const float& b, const float& c, const float& d, std::vector<float> &v)const;
    void quadratic_equation(const float& a, const float& b, const float& c, std::vector<float>& v)const;
    void linear_equation(const float& a, const float& b, std::vector<float>& v)const;
    int sign(const long double& ld)const;
    float landau2D(const float& x, const float& y)const;
    
    
    //Utility Methods
    float sqr(const float& x)const;
    void swap(float& realone, float& realtwo)const;
    
    int nSol_;
    float coeffs_[5];
    std::vector<float> vect_pxv_;
    std::vector<TopSolution> ttSol_;
    TLorentzVector al_;
    TLorentzVector l_;
    TLorentzVector b_;
    TLorentzVector bbar_;
    TLorentzVector top_;
    TLorentzVector topbar_;
    TLorentzVector neutrino_;
    TLorentzVector neutrinobar_;
    TLorentzVector w_;
    TLorentzVector wbar_;
    TLorentzVector tt_;
    
    TLorentzVector true_top_;
    TLorentzVector true_topbar_;
    TLorentzVector true_neutrino_;
    TLorentzVector true_neutrinobar_;
    
    float px_miss_;
    float py_miss_;
    
    float mt_;
    float mtbar_;
    float mb_;
    float mbbar_;
    float mw_;
    float mwbar_;
    float ml_;
    float mal_;
    float mv_;
    float mav_;
    
    float a1_,a2_,a3_,a4_;
    float b1_,b2_,b3_,b4_;
    float c22_,c21_,c20_,c11_,c10_,c00_;
    float d22_,d21_,d20_,d11_,d10_,d00_;
    float d0_,d1_,d2_;
    float c0_,c1_,c2_;
    
};


#endif
