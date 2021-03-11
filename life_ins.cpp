// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include <vector>
#include <boost/optional.hpp>

#define print(var) Rcpp::Rcout<<#var"= "<<var<<"\n"
#define repi(i,a,b) for(int i=int(a);i<=int(b);++i)
#define rrepi(i,a,b) for(int i=int(a);i>=int(b);--i)

using type=double;
using int16=std::int16_t;
type initial = pow(10,11);

using namespace Rcpp;

class base_number
{
private:
public:
  type i,v,q_wi,q_wa;
  int16 x_0,x,m,n,omega;
  std::vector<int16> age;
  std::vector<type> q,l_wi,d_i,d_wi,C_wi,D_wi,M_wi,N_wi
    ,l_wa,d_a,d_wa,C_wa,D_wa,M_wa,N_wa;

  void set_base_rate(const int16 x_0_
                       ,const type i_
                       ,const type q_wi_
                       ,const type q_wa_
                       ,const NumericVector& q_){
    this->x_0=x_0_;
    this->i=i_;
    this->q_wi=q_wi_;
    this->q_wa=q_wa_;
    this->q=as<std::vector<type>>(q_);
  }
  
  void set_policy_data(const int16 x_
                         ,const int16 m_
                         ,const int16 n_){
    this->x=x_;
    this->m=m_;
    this->n=n_;
  }
  
  void calc_base_number(){
    v=1/(1+i);
    omega=x_0+q.size();
    age.resize(omega);
    l_wi.resize(omega);
    d_i.resize(omega);
    d_wi.resize(omega);
    C_wi.resize(omega);
    D_wi.resize(omega);
    M_wi.resize(omega);
    N_wi.resize(omega);
    
    l_wa.resize(omega);
    d_a.resize(omega);
    d_wa.resize(omega);
    C_wa.resize(omega);
    D_wa.resize(omega);
    M_wa.resize(omega);
    N_wa.resize(omega);
    
    repi(t,0,x_0-1){
      age[t]=t;
      l_wi[t]=d_i[t]=d_wi[t]=C_wi[t]=D_wi[t]=M_wi[t]=N_wi[t]=0.;
      l_wa[t]=d_a[t]=d_wa[t]=C_wa[t]=D_wa[t]=M_wa[t]=N_wa[t]=0.;
    }
    repi(t,x_0,omega){
      age[t]=t;
      
      l_wi[t]=(t==x_0?initial:std::max(l_wi[t-1]-d_i[t-1]-d_wi[t-1],0.));
      d_i[t]=l_wi[t]*q[t-x_0]*(1.-1/2.*q_wi/(1.-1/2.*q[t-x_0]));
      d_wi[t]=l_wi[t]*q_wi;
      C_wi[t]=pow(v,t+1/2.)*d_i[t];
      D_wi[t]=pow(v,t)*l_wi[t];
      
      l_wa[t]=(t==x_0?initial:std::max(l_wa[t-1]-d_a[t-1]-d_wa[t-1],0.));
      d_a[t]=l_wa[t]*q[t-x_0]*(1.-1/2.*q_wa/(1.-1/2.*q[t-x_0]));
      d_wa[t]=l_wa[t]*q_wa;
      C_wa[t]=pow(v,t+1/2.)*d_a[t];
      D_wa[t]=pow(v,t)*l_wa[t];
    }
    rrepi(t,omega,x_0){
      M_wi[t]=C_wi[t]+(t==omega?0.:M_wi[t+1]);
      N_wi[t]=D_wi[t]+(t==omega?0.:N_wi[t+1]);
      
      M_wa[t]=C_wa[t]+(t==omega?0.:M_wa[t+1]);
      N_wa[t]=D_wa[t]+(t==omega?0.:N_wa[t+1]);
    }
  }
  
  type a_12(int16 t){
    return (1-pow(v,t))/(1-v);
  }
  
  type a_wi(int16 t){
    if(t<m){
      return (N_wi[x+t]-N_wi[x+m])/D_wi[x+t];
    }else{
      return 0.;
    }
  }
  
  type a_wi_12(int16 t){
    if(t<m){
      return (N_wi[x+t]-N_wi[x+m])/D_wi[x+t]-11/24.*(1-D_wi[x+m]/D_wi[x+t]);
    }else{
      return 0.;
    }
  }
  
  type a_w(int16 t){
    if(t<m){
      return (N_wi[x+t]-N_wi[x+m])/D_wi[x+t]+D_wi[x+m]/D_wi[x+t]*N_wa[x+m]/D_wa[x+m];
    }else if(m<=t){
      return N_wa[x+t]/D_wa[x+t];
    }
  }
  
  type m_a_w(int16 t){
    if(t<m){
      return D_wi[x+m]/D_wi[x+t]*N_wa[x+m]/D_wa[x+m];
    }else if(m<=t){
      return 0.;
    }
  }
  
  DataFrame output_base_number(){
    DataFrame df=DataFrame::create(Named("age")=age
                                     ,Named("l_wi",l_wi)
                                     ,Named("d_i",d_i)
                                     ,Named("d_wi",d_wi)
                                     ,Named("C_wi",C_wi)
                                     ,Named("D_wi",D_wi)
                                     ,Named("M_wi",M_wi)
                                     ,Named("N_wi",N_wi)
                                     ,Named("l_wa",l_wa)
                                     ,Named("d_a",d_a)
                                     ,Named("d_wa",d_wa)
                                     ,Named("C_wa",C_wa)
                                     ,Named("D_wa",D_wa)
                                     ,Named("M_wa",M_wa)
                                     ,Named("N_wa",N_wa)
    );
    return df;
  }
};

base_number p;

// [[Rcpp::export]]
DataFrame test(int x_0,type i,type q_wi,type q_wa,NumericVector q){
  p.set_base_rate(x_0,i,q_wi,q_wa,q);
  p.calc_base_number();
  return p.output_base_number();
}

