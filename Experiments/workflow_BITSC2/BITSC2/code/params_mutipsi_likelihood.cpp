#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
using namespace Rcpp;
using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]



/*
██   ██ ███████  █████  ██████  ███████ ██████
██   ██ ██      ██   ██ ██   ██ ██      ██   ██
███████ █████   ███████ ██   ██ █████   ██████
██   ██ ██      ██   ██ ██   ██ ██      ██   ██
██   ██ ███████ ██   ██ ██████  ███████ ██   ██
*/

arma::mat log_prob_X(arma::mat X, arma::mat D, arma::mat L, arma::mat Z, double w, double mu);
arma::mat log_prob_D(arma::mat D, arma::mat L, arma::vec psi, double rho, double s, double mu);
double prior_C(int k, arma::vec phi);


/* ------------------------ basic ----------------------------------*/
inline bool any_cpp(LogicalVector lv)
{
  return is_true(any(lv));
}

arma::mat matrix_sub(arma::mat M, LogicalVector a, int b)
{
  // b=1: select row
  // b=2: select column
  arma::mat out;
  if(b==2){
    arma::colvec z=as<arma::colvec>(a);
    out=M.cols(find(z==1));
  } else if(b==1){
    arma::rowvec z=as<arma::rowvec>(a);
    out=M.rows(find(z==1));
  }
  return out;
}


/*
████████ ██████  ███████ ███████
   ██    ██   ██ ██      ██
   ██    ██████  █████   █████
   ██    ██   ██ ██      ██
   ██    ██   ██ ███████ ███████
*/


//----------------------------------------------------------------------------------------
/*  -------------------------------- tree related  -------------------------------*/
//----------------------------------------------------------------------------------------

//[[Rcpp::export]]
LogicalVector find_desc_cpp(int k, NumericVector tree)
{
  // returning 1,0 vector indicating whether each subclone is descendant
  LogicalVector parent_set(tree.size());

  parent_set[k-1]=1;
  if (k<tree.size()){
    for (int i=k+1;i<=tree.size();i++){
      parent_set[i-1] = (parent_set[tree[i-1]-1]) ? 1:0;
    }
  }
  return parent_set;
}

//[[Rcpp::export]]
LogicalVector find_ances_cpp(int k, NumericVector tree)
{
  // returning 1,0 vector indicating whether each subclone is ancestor
  LogicalVector child_set(tree.size());

  child_set[k-1]=1;
  int kk=k;
  do{
    child_set[tree[kk-1]-1]=1;
    kk=tree[kk-1];
  }while(kk!=0);
    
  return child_set;
}


//[[Rcpp::export]]
LogicalVector find_child_cpp(double k, NumericVector tree)
{
  return tree==k;
}


//[[Rcpp::export]]
NumericMatrix Lo_to_L_cpp(NumericMatrix B, NumericVector tree)
{
  int M=B.nrow();// number of muts
  int K=tree.length(); // number of subclones

  NumericMatrix out(M,K);
  NumericVector a(2);
  NumericVector b(K);

  LogicalVector temp(K);
  for (int i=0;i<M;i++){
    a=B.row(i);
    b=wrap(2*arma::ones(K));

    if(a[0]!=0){
      temp=find_desc_cpp(a[0], tree);
      b[temp] = a[1]+2;
    }

    out.row(i)=b;
  }
  return out;
}

//[[Rcpp::export]]
NumericMatrix L_to_Lo_cpp(NumericMatrix B)
{
  int M=B.nrow();// number of muts
  int K=B.ncol(); // number of subclones

  NumericMatrix out(M,2);

  NumericVector temp(2);
  for (int m=0;m<M;m++){
    temp[0]=temp[1]=0;
      for(int k=0;k<K;k++){
        if(B(m,k) != 2){temp[1]= B(m,k)-2;temp[0]=k+1;break;}
      }
    out.row(m)=temp;
  }
  return out;
}


//[[Rcpp::export]]
NumericMatrix Zo_to_Z_cpp(NumericMatrix Zo, NumericMatrix Lo, NumericVector g, NumericVector tree)
{
  int M=Zo.nrow();// number of muts
  int K=tree.length(); // number of subclones

  NumericMatrix out(M,K);
  NumericVector a(2);
  NumericVector b(K);
  LogicalVector temp1(K);
  LogicalVector temp2(K);

  for (int i=0;i<M;i++){
    a=Zo.row(i);
    b=wrap(arma::zeros(K));
    temp1=find_desc_cpp(a[0], tree);
    b[temp1] = a[1];
    if(g.at(i)){
      temp2=find_desc_cpp(Lo.at(i,0), tree);
      b[temp2] = a[1] + Lo.at(i,1);
    }
    out.row(i)=b;
  }
  return out;
}

//[[Rcpp::export]]
NumericMatrix Z_to_Zo_cpp(NumericMatrix Z)
{
  int M=Z.nrow();// number of muts
  int K=Z.ncol(); // number of subclones

  NumericMatrix out(M,2);

  NumericVector temp(2);
  for (int m=0;m<M;m++){
    temp[0]=temp[1]=0;
    for(int k=0;k<K;k++){
      if(Z(m,k) != 0){temp[1]= 1;temp[0]=k+1;break;}
    }
    out.row(m)=temp;
  }
  return out;
}



/*
███████  █████  ███    ███ ██████  ██      ███████      ██████  ██   ██ ████████
██      ██   ██ ████  ████ ██   ██ ██      ██           ██   ██ ██   ██    ██
███████ ███████ ██ ████ ██ ██████  ██      █████        ██████  ███████    ██
     ██ ██   ██ ██  ██  ██ ██      ██      ██           ██      ██   ██    ██
███████ ██   ██ ██      ██ ██      ███████ ███████      ██      ██   ██ ████████
*/

//------------------------------------------------------------------------
// ------------------------- sampling of theta --------------------------
//------------------------------------------------------------------------


//[[Rcpp::export]]
double prior_theta(arma::vec theta, List Params)
{
  double r=Params["r"];
  double out= (r-1)*sum(log(theta)) - sum(theta);
  return out;
}

//[[Rcpp::export]]
double log_p_theta_cpp(arma::vec theta, List samp, List Params, double temper)
{
  double out=0;
  int N=Params["N"];
  NumericVector C=samp["C"];
  double G = sum(theta);
  arma::vec phi = theta/G;
  
  for(int i=0;i<N;i++){
    int k = C.at(i);
    out += prior_C(k,phi);
  }
  out += prior_theta(theta, Params);
  return out/temper;
}


NumericVector samp_theta_cpp(arma::vec theta, int i, List samp, List Params, double p0, double tune_par, double temper, double err)
{
  NumericVector out(2);

  double theta0= theta[i];

  //propose
  double theta1= R::rgamma(tune_par*theta0, 1/tune_par);
  if(theta1==0){theta1 += err;}



  double p_01= R::dgamma(theta1,tune_par*theta0,1/tune_par,0);
  double p_10= R::dgamma(theta0,tune_par*theta1,1/tune_par,0);


  //posterior of new data

  if(p0==0){p0=log_p_theta_cpp(theta,samp,Params,temper);}
  arma::vec new_theta = theta;
  new_theta[i]=theta1;

  double p1=log_p_theta_cpp(new_theta,samp,Params,temper);
  // acceptance probability

  double temp_prob=exp(p1-p0)*(p_10/p_01);
  double acc_prob= (temp_prob >= 1)? 1 : temp_prob;


  double u=R::runif(0,1);
  if(u<=acc_prob){
    out[0]=theta1;
    out[1]=p1;
  } else{
    out[0]=theta0;
    out[1]=p0;
  }

  return out;
}

//[[Rcpp::export]]
arma::vec samp_theta_all (List samp, List Params, double theta_tune,double temper)
{
  arma::vec out= samp["theta"];
  int K=out.n_rows;
  double p0=0;

  for (int i=0; i<K;i++){
    NumericVector update= samp_theta_cpp(out,i,samp,Params,p0,theta_tune,temper,0.01);
    out.at(i)=update[0];
    p0=update[1];
    // update vector
    if(i==(K-1)) {p0=0;}
  }
  return out;
}


/*
███████  █████  ███    ███ ██████  ██      ███████      ███████
██      ██   ██ ████  ████ ██   ██ ██      ██          ██
███████ ███████ ██ ████ ██ ██████  ██      █████       ██
     ██ ██   ██ ██  ██  ██ ██      ██      ██          ██
███████ ██   ██ ██      ██ ██      ███████ ███████      ███████
*/



//------------------------------------------------------------------------
// ------------------------- sampling of C -------------------------------
//------------------------------------------------------------------------

//[[Rcpp::export]]
double prior_C(int k, arma::vec phi)
{
  double out = phi.at(k-1);
  return log(out);
}


//[[Rcpp::export]]
double log_p_C_cpp(int n, int k, arma::vec Zk, arma::vec Lk, arma::mat X, arma::mat D, List Params, List samp, double temper)
{
  double out;
  arma::vec psi = as<arma::vec>(Params["psi"]);
  arma::vec phi = samp["phi"];
  double rho = samp["rho"];
  double w = samp["w"];
  double mu = samp["mu"];
  double s = samp["s"];
  
  int M = Params["M"];
  arma::vec psi_n(M);
  psi_n.fill(psi[n-1]);
  
  arma::mat px = log_prob_X(X.col(n-1), D.col(n-1), Lk, Zk, w, mu);
  arma::mat pd = log_prob_D(D.col(n-1), Lk, psi_n, rho, s,mu);
  out =  accu(px) + accu(pd) + prior_C(k, phi);
  return out/temper;
}

//[[Rcpp::export]]
int Gibbs_C_cpp(int n, arma::mat Z, arma::mat L, arma::mat X, arma::mat D, List Params, List samp, double temper){
  int K=Params["K"];
  arma::vec p_vec=arma::zeros<arma::vec>(K);
  arma::vec Zk,Lk;
  for(int k=1; k<=K; k++){
    Zk=Z.col(k-1);
    Lk=L.col(k-1);
    p_vec.at(k-1)= log_p_C_cpp(n, k, Zk, Lk, X, D, Params, samp, temper);
  }
  double increase = max(p_vec);
  int s= numel(p_vec);
  arma::vec p_vec_max(s);
  p_vec_max.fill(increase);
  arma::vec p_vec_norm = exp(p_vec - p_vec_max);
  
  IntegerVector frame=Range(1,K);
  int out= as<int>(Rcpp::RcppArmadillo::sample(frame,1,1,p_vec_norm));
  return out;
}


//[[Rcpp::export]]
NumericVector samp_C_all(List samp, List Params, arma::mat X, arma::mat D, double temper)
{
  NumericVector out=samp["C"];
  int N=Params["N"], K=Params["K"];
  arma::mat Z=samp["Z"], L=samp["L"];
  for (int n=1; n<=N; n++){
    int temp=Gibbs_C_cpp(n,Z,L,X,D,Params,samp,temper);
    out(n-1)=temp;
  }
  return out;
}

/*
███████  █████  ███    ███ ██████  ██      ███████     ██
██      ██   ██ ████  ████ ██   ██ ██      ██          ██
███████ ███████ ██ ████ ██ ██████  ██      █████       ██
     ██ ██   ██ ██  ██  ██ ██      ██      ██          ██
███████ ██   ██ ██      ██ ██      ███████ ███████     ███████
*/



//----------------------------------------------------------------
// ------------------------- sampling of L--------------------------
//--------------------------------------------------------------------

//[[Rcpp::export]]
double prior_L_cpp(NumericVector Lm, List samp, List Params)
{
  double out;

  double K=Params["K"];
  double pi=samp["pi"];
  double max_CN=Params["max_CN"];

  LogicalVector temp= (Lm != 2);
  if(any_cpp(temp)){
    out=(1-pi)/((K-1)*max_CN);
  } else {
    out=pi;
  }
  return log(out);
}


//[[Rcpp::export]]
double log_p_Lo_cpp(arma::vec Lm,arma::vec Zm, int m, arma::mat D, arma::mat X, List samp, List Params, double temper)
{
  int N = Params["N"];
  arma::vec psi = as<arma::vec>(Params["psi"]);
  NumericVector C = samp["C"];
  double rho = samp["rho"];
  double w = samp["w"];
  double mu = samp["mu"];
  double s = samp["s"];

  arma::vec Ln(N);
  arma::vec Zn(N);
  for(int i=0;i<N;i++){
    Ln.at(i)=Lm.at(C.at(i)-1);
    Zn.at(i)=Zm.at(C.at(i)-1);
  }
  arma::mat px = log_prob_X(X.row(m-1), D.row(m-1), Ln.t(), Zn.t(), w, mu);
  arma::mat pd = log_prob_D(D.row(m-1), Ln.t(), psi, rho, s, mu);
  
  double out = accu(pd) + accu(px);

  out /= temper;
  return out;
}


//[[Rcpp::export]]
NumericVector propose_Lo_cpp(NumericVector tree, arma::mat Z_seg, List Params)
{
  int K=Params["K"];
  int max_CN=Params["max_CN"];

  int origin=R::runif(2,K+1); // randomly select an origin
  LogicalVector desc=find_desc_cpp(origin, tree);

  arma::mat temp=matrix_sub(Z_seg,desc,2);
  int min_cnv=temp.max(); // make sure L>Z
  int new_cn= R::runif(min_cnv,max_CN+1); // new copy number

  NumericVector out(2);
  out.at(0)=origin;
  out.at(1)=new_cn-2;

  return out;
}



//[[Rcpp::export]]
NumericVector Gibbs_Lo_cpp(arma::uvec indVec, List samp, arma::mat X, arma::mat D, List Params, double temper)
{
  arma::mat L=samp["L"], Z=samp["Z"], Zo=samp["Zo"];
  int K= Params["K"];
  int maxCN=Params["max_CN"];
  NumericVector tree=samp["Ttree"];
  arma::mat Z_seg = Z.rows(indVec-1);
  NumericVector g(1);

  arma::vec p_vec=arma::zeros<arma::vec>((K-1)*maxCN+1);
  double p0;

  // k=0
  int min_cnv=Z_seg.max(); // make sure L>Z
  arma::vec normal_l(K);
  arma::vec normal_lo(2);
  normal_l.fill(2);
  normal_lo.at(0)=0;
  normal_lo.at(1)=0;

  if(min_cnv<=2){
    p0= prior_L_cpp(wrap(normal_l),samp,Params)/temper;
    for(arma::uvec::iterator it_m=indVec.begin(); it_m != indVec.end(); it_m++){
      arma::vec Zm = Z.row(*it_m-1).t();
      p0 += log_p_Lo_cpp(normal_l,Zm,*it_m,D,X,samp,Params,temper);
    }
    p_vec.at(0)=p0;
  }


  int ind=1;
  arma::vec new_l(K);
  for (int k=2; k<=K; k++){
    for(int v=-2; v<=maxCN-2 ; v++){
      if(v==0) continue;
      LogicalVector desc=find_desc_cpp(k, tree);
      LogicalVector cnv_ances=find_ances_cpp(k,tree);
      arma::vec cur_lo(2);
      cur_lo.at(0)= k;
      cur_lo.at(1)= v;
      
      min_cnv=matrix_sub(Z_seg,desc,2).max();
      if(v+2>=min_cnv){
         // create l vector
         new_l.fill(2);
         new_l.elem(find(as<arma::uvec>(desc))) = new_l.elem(find(as<arma::uvec>(desc))) + v;
         // calculate posterior
         p0=prior_L_cpp(wrap(new_l),samp,Params)/temper;
        
         for(arma::uvec::iterator it_m=indVec.begin(); it_m != indVec.end(); it_m++){
           arma::vec Zm = Z.row(*it_m-1).t();
           NumericMatrix Zmo = Z_to_Zo_cpp(wrap(Zm.t()));
           if(-v>Zmo.at(0,1)&&-v>2+v-Zmo.at(0,1)){
             p0=0;
             break;
           }
           
           g.at(0)=0;
           NumericMatrix temp_Z = Zo_to_Z_cpp(Zmo,wrap(cur_lo.t()),g,tree);
           double p1 = log_p_Lo_cpp(new_l,as<arma::vec>(temp_Z),*it_m,D,X,samp,Params,temper);
           
           if(cnv_ances[Zmo.at(0,0)-1]){
             NumericMatrix temp_Z1;
             double p2;
          
             if(v>0){
               g.at(0)=1;
               temp_Z1 = Zo_to_Z_cpp(Zmo,wrap(cur_lo.t()),g,tree);
               p2 = log_p_Lo_cpp(new_l,as<arma::vec>(temp_Z1),*it_m,D,X,samp,Params,temper);
             }
             if(v<0){
               if(-v<=Zmo.at(0,1)){
                 g.at(0)=1;
                 temp_Z1 = Zo_to_Z_cpp(Zmo,wrap(cur_lo.t()),g,tree);
                 p2 = log_p_Lo_cpp(new_l,as<arma::vec>(temp_Z1),*it_m,D,X,samp,Params,temper);
               }else{
                 p2=p1;
               }
             }

             p0 = p0 + (log(exp(p1*temper)+exp(p2*temper))-log(2))/temper;
           }
           if(!cnv_ances[Zmo.at(0,0)-1]){
             p0=p0+p1;
           }
         }
        if(p0<(-pow(10,100))){
          p0 = -pow(10,100);
        }
         p_vec.at(ind)=p0;
       }
       ind++;
    }
  }
  

  double increase= p_vec.elem(find(p_vec)).max();
  p_vec.elem(find(p_vec))= exp( p_vec.elem(find(p_vec)) - increase);
  
  IntegerVector frame=Range(0,(K-1)*maxCN);
  int ind_samp= as<int>(Rcpp::RcppArmadillo::sample(frame,1,1,p_vec));
  if (ind_samp==0) {return wrap(normal_lo);}
  

  int k=2+(ind_samp-1)/maxCN;
  int v=(ind_samp-1)-(k-2)*maxCN-2;
  v= (v>=0) ? v+1 : v;

  arma::vec new_lo(2);
  new_lo.at(0)= k;
  new_lo.at(1)= v;

  return wrap(new_lo);
}

//[[Rcpp::export]]
NumericMatrix samp_Lo_all(List samp, arma::mat X, arma::mat D, List Params, double temper)
{
  int M=Params["M"]; // number locus
  int K=Params["K"]; // number of subclones
  arma::mat segs=Params["segments"];
  int nsegs= segs.n_rows;

  NumericMatrix out(M,2);
  NumericVector temp(2);

  for (int i=0; i<nsegs; i++){
    // generate sequence for this segment
    int start=segs(i,0), end=segs(i,1);
    arma::uvec one_seg(end-start+1);
    for(arma::uvec::iterator uvec_it=one_seg.begin(); uvec_it != one_seg.end(); uvec_it++){
      *uvec_it=start++;
    }
    temp=Gibbs_Lo_cpp(one_seg,samp,X,D,Params,temper);

    // assign output to all rows in this segment
    for(arma::uvec::iterator uvec_it=one_seg.begin(); uvec_it != one_seg.end(); uvec_it++){
      out.row(*uvec_it-1)=temp;
    }
  }

  return out;
}


/*
███████  █████  ███    ███ ██████  ██      ███████     ███████
██      ██   ██ ████  ████ ██   ██ ██      ██             ███
███████ ███████ ██ ████ ██ ██████  ██      █████         ███
     ██ ██   ██ ██  ██  ██ ██      ██      ██           ███
███████ ██   ██ ██      ██ ██      ███████ ███████     ███████
*/



//---------------------------------------------------------------------------
// ------------------------- sampling of Z  --------------------------
//------------------------------------------------------------------------
//[[Rcpp::export]]
double prior_Z_cpp(NumericVector Zm, List Params)
{
  int a=max(Zm);
  int max_mut=Params["max_mut"];
  double zeta=Params["zeta"];

  double p_a=pow(zeta,a);
  double p_all=(1-pow(zeta,max_mut+1))/(1-zeta);
  p_all -= 1;

  double out=p_a/p_all;
  return log(out);
}


//[[Rcpp::export]]
double log_p_Z_cpp(int m, arma::vec Zm, arma::vec Lm,  arma::mat D, arma::mat X, List samp, List Params, double temper)
{
  int N = Params["N"];
  NumericVector C = samp["C"];
  double w = samp["w"];
  double mu = samp["mu"];
  
  arma::vec Ln(N);
  arma::vec Zn(N);
  for(int i=0;i<N;i++){
    Ln.at(i) = Lm.at(C.at(i)-1);
    Zn.at(i) = Zm.at(C.at(i)-1);
  }
  
  arma::mat px = log_prob_X(X.row(m-1),D.row(m-1),Ln.t(),Zn.t(),w,mu);
  double out = accu(px) + prior_Z_cpp(wrap(Zm),Params);

  out = out/temper;
  return out;
}



//[[Rcpp::export]]
NumericVector propose_Zo_cpp(NumericVector tree, NumericVector Lm, List Params)
{
  int K=Params["K"];
  // sample an origin with L > 0
  LogicalVector loc=(Lm > 0);
  loc[0]=0; // change value for normal clone

  int s=sum(loc);
  int u=R::runif(1,s+1);
  int cumsum=0, index=1;
  for (LogicalVector::iterator it=loc.begin(); it !=loc.end(); it++){
    cumsum += loc[index-1];
    if(cumsum==u)break;
    index++;
  }

  LogicalVector desc=find_desc_cpp(index,tree); // subclones to be changed

  //sample number of mutant copy
  NumericVector Lm_desc=Lm[desc];
  int min_Lm_desc= min(Lm_desc);
  NumericVector out(2);
  
  if(min_Lm_desc>0){
    //int max_mut2 = (max_mut > min_Lm_desc) ? min_Lm_desc : max_mut; // maximum mutant copy allowed
    //u=R::runif(1,max_mut2+1);
    out.at(0)=index;
    out.at(1)=1;
  }
  return out;
}



//[[Rcpp::export]]
NumericVector Gibbs_Zo_cpp(int m, arma::vec Lmo, List samp, arma::mat X, arma::mat D, List Params, double temper)
{
  NumericVector tree = samp["Ttree"];
  int K = Params["K"];
  int maxMut = Params["max_mut"];
  
  NumericVector Lm1 = Lo_to_L_cpp(wrap(Lmo.t()),tree);
  arma::vec Lm = as<arma::vec>(Lm1);
  
  LogicalVector cnv_ances = find_ances_cpp(Lmo.at(0),tree);
  arma::vec p_vec = arma::zeros<arma::vec>((K-1)*maxMut);
  int ind=0;
  NumericVector g(1);
  
  for(int k=2; k<=K; k++){
    for(int v=1; v<=maxMut; v++){
      //prob for (k,v)
      arma::uvec q=find(as<arma::uvec>(find_desc_cpp(k,tree)));
      int temp_max=Lm.elem(q).min();
      
      if(temp_max<v) {ind++;continue;}
      
      arma::vec cur_zo(2);
      cur_zo.at(0) = k;
      cur_zo.at(1) = v;
      
      
      g.at(0)=0;
      double p1;
      NumericMatrix temp_Z = Zo_to_Z_cpp(wrap(cur_zo.t()),wrap(Lmo.t()),g,tree);
      p1 = log_p_Z_cpp(m,as<arma::vec>(temp_Z),Lm,D,X,samp,Params,temper);
      
      
      if(cnv_ances[k-1]){
        double p2;
        if(Lmo.at(1)>0){
          g.at(0)=1;
          NumericMatrix temp_Z = Zo_to_Z_cpp(wrap(cur_zo.t()),wrap(Lmo.t()),g,tree);
          p2 = log_p_Z_cpp(m,as<arma::vec>(temp_Z),Lm,D,X,samp,Params,temper);
        }
        if(Lmo.at(1)<0){
          if(v>=-Lmo.at(1)){
            g.at(0)=1;
            NumericMatrix temp_Z = Zo_to_Z_cpp(wrap(cur_zo.t()),wrap(Lmo.t()),g,tree);
            p2 = log_p_Z_cpp(m,as<arma::vec>(temp_Z),Lm,D,X,samp,Params,temper);
          }else{
            p2 = p1;
          }
        }

        double pp=(log(exp(p1*temper)+exp(p2*temper))-log(2))/temper;
        if(pp<(-pow(10,100))){
          pp = -pow(10,100);
        }
        p_vec.at(ind)=pp;
      }
      
      if(!cnv_ances[k-1]){
        p_vec.at(ind) = p1;
      }
      ind++;
    }
  }
  
  double increase = p_vec.elem(find(p_vec)).max();
  p_vec.elem(find(p_vec)) = exp( p_vec.elem(find(p_vec)) - increase);
  IntegerVector frame = Range(1,(K-1)*maxMut);
  int ind_samp= as<int>(Rcpp::RcppArmadillo::sample(frame,1,1,p_vec));
  
  arma::vec new_zo(2);
  int k = 2+(ind_samp-1)/maxMut;
  int v = ind_samp-(k-2)*maxMut;
  
  new_zo.at(0)= k;
  new_zo.at(1)= v;
  
  return wrap(new_zo);
}



//[[Rcpp::export]]
NumericMatrix samp_Zo_all(List samp, arma::mat X, arma::mat D, List Params, double temper)
{
  int M=Params["M"], K=Params["K"];
  arma::mat Lo=samp["Lo"];
  NumericMatrix out(M,2);
  NumericVector temp(2);

  arma::vec Zmo,Lmo;
  for (int m=1; m<=M; m++){
    Lmo = Lo.row(m-1).t();
    temp = Gibbs_Zo_cpp(m,Lmo,samp,X,D,Params,temper);
    out.row(m-1)=temp;
  }
  return out;
}


//[[Rcpp::export]]
NumericVector estimate_g(List samp, arma::mat X, arma::mat D, List Params)
{
  int N = Params["N"];
  int M = Params["M"];
  double w = samp["w"];
  double mu = samp["mu"];
  NumericVector tree = samp["Ttree"];
  arma::mat L = samp["L"], Zo = samp["Zo"];
  NumericVector C = samp["C"];
  NumericMatrix Lo = samp["Lo"];

  arma::vec Zmo,Lmo;
  NumericVector Lm;
  NumericVector g(M);
  NumericMatrix Zm0, Zm1;
  NumericVector gm(1);
  
  for(int m=1; m<=M; m++){
    Lm = L.row(m-1).t();
    Lmo = Lo.row(m-1);
    Zmo = Zo.row(m-1).t();
    
    LogicalVector cnv_ances = find_ances_cpp(Lmo.at(0),tree);
    
    if(cnv_ances[Zmo.at(0)-1]){
      gm.at(0)=0;
      Zm0 = Zo_to_Z_cpp(wrap(Zmo.t()),wrap(Lmo.t()),gm,tree);
      gm.at(0)=1;
      Zm1 = Zo_to_Z_cpp(wrap(Zmo.t()),wrap(Lmo.t()),gm,tree);
      arma::vec Ln(N);
      arma::vec Zn0(N),Zn1(N);
      for(int i=0;i<N;i++){
        Ln.at(i) = Lm.at(C.at(i)-1);
        Zn0.at(i) = Zm0.at(0,C.at(i)-1);
        Zn1.at(i) = Zm1.at(0,C.at(i)-1);
      }
      
      double p0 = accu(log_prob_X(X.row(m-1),D.row(m-1),Ln.t(),Zn0.t(),w,mu));
      double p1 = accu(log_prob_X(X.row(m-1),D.row(m-1),Ln.t(),Zn1.t(),w,mu));

      if(p0>=p1){
        g.at(m-1)=0;
      }
      if(p0<p1){
        g.at(m-1)=1;
      }
    }else{
      g.at(m-1)=0;
    }
  }
  return g;
}

/*
 ██████  ████████ ██   ██ ███████ ██████  ███████
██    ██    ██    ██   ██ ██      ██   ██ ██
██    ██    ██    ███████ █████   ██████  ███████
██    ██    ██    ██   ██ ██      ██   ██      ██
 ██████     ██    ██   ██ ███████ ██   ██ ███████
*/


// ------------------------------------------------------------
/*  --------------------- others ---------------------- */
// ------------------------------------------------------------

// [[Rcpp::export]]
double corr_pro(double p)
{
  double epsilon = pow(10,-3);
  if(p>0 && p<1){
    p = p-2*epsilon/3;
  }
  if(p==1){
    p=1-epsilon;
  }
  if(p==0){
    p=epsilon;
  }
  return p;
}

// [[Rcpp::export]]
double dbbC( int x,int N, double u, double v, bool log)
{
  double logval = R::lbeta(x + u, N - x + v) - R::lbeta(u, v) + R::lchoose(N, x);
  double ret;
    if (log) {
        ret = logval;
    }
    else {
        ret = exp(logval);
    }
    return(ret);
}



//[[Rcpp::export]]
arma::mat log_prob_X(arma::mat X, arma::mat D, arma::mat L, arma::mat Z, double w, double mu)
{
  int M = X.n_rows;
  int N = X.n_cols;
  double epsilon = pow(10,-3);
  double u;
  double v;
  double u1 = epsilon*w;
  double v1 = (1-epsilon)*w;
  // arma::mat P = Z/L;
  arma::mat out(M,N);
  
  for (int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      
      if(L.at(i,j)==0){
        if(Z.at(i,j)==0){
          double b1 = dbbC(X.at(i,j),D.at(i,j),u1,v1,0);
          out.at(i,j) = b1;
        }
        else{
          out.at(i,j) = 0;
        }
      }
      else{
        double p = Z.at(i,j)/L.at(i,j);
        if(p>1){
          out.at(i,j) = 0;
        }
        else{
          p = corr_pro(p);
          u = p*w;
          v = (1-p)*w;
          
          if(L.at(i,j)==1){
            double b1 = dbbC(X.at(i,j),D.at(i,j),u1,v1,0);
            double b = dbbC(X.at(i,j),D.at(i,j),u,v,0);
            out.at(i,j) = mu * b1 + (1-mu) * b;
          }
          if(L.at(i,j)>1){
            if(Z.at(i,j)==0){
              double b1 = dbbC(X.at(i,j),D.at(i,j),u1,v1,0);
              out.at(i,j) = b1;
            }
            if(Z.at(i,j)>0){
              double p2 = (Z.at(i,j)-1)/(L.at(i,j)-1);
              p2 = corr_pro(p2);
              
              double u2 = p2*w;
              double v2 = (1-p2)*w;
              
              double b2 = dbbC(X.at(i,j),D.at(i,j),u2,v2,0);
              double b = dbbC(X.at(i,j),D.at(i,j),u,v,0);
              
              if(Z.at(i,j)==L.at(i,j)){
                out.at(i,j) = mu* b2 + (1-mu) * b;
              }else{
                double p3 = Z.at(i,j)/(L.at(i,j)-1);
                p3 = corr_pro(p3);
                double u3 = p3*w;
                double v3 = (1-p3)*w;
                
                double b3 = dbbC(X.at(i,j),D.at(i,j),u3,v3,0);
                
                out.at(i,j) = mu*(Z.at(i,j)/L.at(i,j) * b2 + (L.at(i,j)-Z.at(i,j))/L.at(i,j) * b3) + (1-mu) * b;
              }
            }
          }
        }
      }
    }
  }
  for (int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      if(out.at(i,j)==0){
        out.at(i,j)=-999999999;
      }else{
        out.at(i,j)=log(out.at(i,j));
      }
    }
  }
  return out;
}

//[[Rcpp::export]]
double delta_D(double d){
  double delta_d;
  if(d==0){
    delta_d=1;
  }else{
    delta_d=0;
  }
  return(delta_d);
}


//[[Rcpp::export]]
arma::mat log_prob_D(arma::mat D, arma::mat L, arma::vec psi, double rho, double s, double mu)
{
  int M = D.n_rows;
  int N = D.n_cols;
  arma::mat out(M,N);
  double epsilon = pow(10,-3);
  double lambda0 = s*epsilon/(1-epsilon);
  for (int i=0; i<M; i++){
    for(int j=0; j<N; j++){
      double lambda = psi[j] * L.at(i,j)/2;
      if(L.at(i,j)==0){
        out.at(i,j) = log(rho * delta_D(D.at(i,j))+(1-rho) * R::dnbinom_mu(D.at(i,j),s,lambda0,0));
      }else{
        double lambda1 = psi[j]*(L.at(i,j)-1)/2;
        if(L.at(i,j)==1){
          out.at(i,j) = log(rho * delta_D(D.at(i,j))+(1-rho) * (mu * R::dnbinom_mu(D.at(i,j),s,lambda0,0) + (1-mu) * R::dnbinom_mu(D.at(i,j),s,lambda,0)));
        }else{
          out.at(i,j)= log(rho * delta_D(D.at(i,j))+(1-rho) * (mu * R::dnbinom_mu(D.at(i,j),s,lambda1,0) + (1-mu) * R::dnbinom_mu(D.at(i,j),s,lambda,0)));
        }
      }
      if(out.at(i,j)<-999999999){
        out.at(i,j) = -999999999;
      }
    }
  }
  return out;
}


// calculate prior for all parameters
//[[Rcpp::export]]
double log_prior_all(List samp, List Params, double temper)
{
  int N=Params["N"];
  int M=Params["M"];
  arma::mat L=samp["L"], Z=samp["Z"];
  arma::vec theta=samp["theta"], phi=samp["phi"];
  NumericVector C=samp["C"];
  double rho = samp["rho"];
  double mu = samp["mu"];
  double w = samp["w"];
  double s = samp["s"];
  // double omega = samp["over_omega"];

  int a=Params["a_pi"];
  int b=Params["b_pi"];
  double pi=samp["pi"];

  double out=0;
// pi
  double temp=R::dbeta(pi,a,b,1);
  out += temp;
  
//theta prior
  out += prior_theta(theta, Params);
  
//w prior
  double ws_rate=Params["ws_rate"];
  int scale = (int)(1/ws_rate);
  out += R::dgamma(w,Params["ws_shape"], scale, 1);

//s prior
  out += R::dgamma(s,Params["ws_shape"], scale, 1);
  
//C prior
  for (int n=0; n<N; n++){
    out += prior_C(C.at(n), phi);
  }

//L prior
  for(int m=0; m< M; m++){
    out += prior_L_cpp(wrap(L.row(m).t()),samp,Params);
  }

//Z prior
  for(int m=0; m<M; m++){
    out += prior_Z_cpp(wrap(Z.row(m).t()),Params);
  }

// omega prior
  // out +=  log_prior_omega(omega, Params);

  return out/temper;
}


// calculate likilihood matrix l_mn=p(x_mn | ...) p(d_mn | ...)
//[[Rcpp::export]]
arma::mat likelihood_mat(arma::mat X, arma::mat D, arma::mat L, arma::mat Z, double w,double s,double mu, List Params, double rho)
{
  arma::vec psi=Params["psi"];
  
  //arma::mat P= Z / L;
  arma::mat out = log_prob_X(X,D,L,Z,w,mu) + log_prob_D(D,L,psi,rho,s,mu);

 return out;
}
