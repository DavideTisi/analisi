#include "cg.h"
#include "config.h"
#include "../rnd.h"
#include "../lib/json.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <functional>

#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#else
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif


#define NATOMS 36
using json = nlohmann::json;

template <int D>
class QuadraticForm : public Function <Eigen::Matrix<double,D,1>, double> {
  public:
    virtual double operator () (const Eigen::Matrix<double,D,1> & x) final {
        return (0.5*x.transpose()*A*x-b.transpose()*x)(0,0);
    }

    virtual Eigen::Matrix<double,D,1> deriv (const Eigen::Matrix<double,D,1> & x) final {
        return (A*x - b).eval();
    }

    void set_A_b(const Eigen::Matrix<double,D,D> & A_, const Eigen::Matrix<double,D,1> &b_) {
        A=A_;
        b=b_;
        x0=A.colPivHouseholderQr().solve(b);
        std::cout<<"sol:" << x0<<"\n";
    }
    double get_solution_distance(const Eigen::Matrix<double,D,1> & x) {
        return (x-x0).norm();
    }
protected:
    Eigen::Matrix<double,D,D> A;
    Eigen::Matrix<double,D,1> b,x0;
};


template <int N,int DIM,int flags >
class MultiPair : public Function <Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> >, double,  Eigen::Matrix<double,N*DIM,N*DIM>,  Eigen::Matrix<double,N*DIM,1>,
       Eigen::Ref<Eigen::Matrix<double,N*DIM,1> > > {
public:
    virtual double operator() (const Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> > & x ) final {
        double res=0.0;
        for (unsigned int i=0;i<N;i++) {
            for (unsigned int j=i+1;j<N;j++) {
                double r2;
                Eigen::Matrix<double,DIM,1> dx;
                if (want_pbc()){
                    dx=pbc(x.template segment<DIM>(i*DIM,DIM),x.template segment<DIM>(j*DIM,DIM),r2);
                }else{

                    dx=(x.template segment<DIM>(i*DIM,DIM)-x.template segment<DIM>(j*DIM,DIM));
                    r2=dx.squaredNorm();
                }
                res+=pair(r2);
            }
        }
        return res;
    }
    virtual Eigen::Matrix<double,N*DIM,1> deriv(const Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> > & x) final {
        Eigen::Matrix<double,N*DIM,1> res=Eigen::Matrix<double,N*DIM,1>::Zero();
        for (unsigned int i=0;i<N;i++) {
            for (unsigned int j=i+1;j<N;j++) {
                if (i==j)
                    continue;
                double r2;
                Eigen::Matrix<double,DIM,1> dx;
                if (want_pbc()){
                    dx=pbc(x.template segment<DIM>(i*DIM,DIM),x.template segment<DIM>(j*DIM,DIM),r2);
                }
                else {
                    dx=(x.template segment<DIM>(i*DIM,DIM)-x.template segment<DIM>(j*DIM,DIM));
                    r2=dx.squaredNorm();
                }
                double dp_dr=pair_deriv_r2(r2);
                for (unsigned int i0=0;i0<DIM;i0++) {
                    res(i*DIM+i0)+=dp_dr*dx(i0);
                    res(j*DIM+i0)-=dp_dr*dx(i0);
                }
            }
        }
        return res;
    }

    virtual Eigen::Matrix<double,N*DIM,N*DIM> hessian_deriv(const Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> > & x, Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > res1 ) final {
        Eigen::Matrix<double,N*DIM,N*DIM> res;
//        Eigen::Matrix<double,N*DIM,N*DIM> res1;
        if (has_deriv2()) {
            for (unsigned int i1=0;i1<N;i1++) {
                for (unsigned int i2=i1+1;i2<N;i2++) { //hessian is symmetric
                    // note: here I calculate also the diagonal term, that has a different form with one more sum
                    // If one makes the calculation,
                    // the diagonal term is simply the sum of the other terms in the row. Because the matrix is symmetric, at the end of the day,
                    // for every out of diagonal term (i1,i2) that I calculate, I add the same quantity to the diagonal terms (i1,i1) and (i2,i2)
                    // the factor in front of the 3x3 matrix is the same because of the symmetry of the potential

                    //get pbc distance
                    Eigen::Matrix<double,DIM,1> dx;
                    double r2;
                    if (want_pbc()){
                        dx=pbc(x.template segment<DIM>(i1*DIM,DIM),x.template segment<DIM>(i2*DIM,DIM),r2);
                    }
                    else {
                        dx=(x.template segment<DIM>(i1*DIM,DIM)-x.template segment<DIM>(i2*DIM,DIM));
                        r2=dx.squaredNorm();
                    }

                    double f2=4*pair_deriv2_r2(r2); //second derivative with respect to r^2_ij
                    double f1=2*pair_deriv_r2(r2); //first derivative with respect to r^2_ij
                    double sub[DIM*DIM]={0.0};
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        // calculate also force (it is free, at this point)
                        res1(i1*DIM+d0)+=f1*dx(d0);
                        res1(i2*DIM+d0)-=f1*dx(d0);

                        // diagonal only term
                        sub[d0*DIM+d0]+=f1;
                        for (unsigned int d1=d0;d1<DIM;d1++) {
                            sub[d0*DIM+d1]+=f2*dx(d0)*dx(d1);
                        }
                    }
                    /*
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        for (unsigned int d1=0;d1<DIM;d1++)
                            std::cout << sub[d0*DIM+d1]<<" ";
                        std::cout << "\n";
                    }
                    std::cout << "\n\n";
*/
                    //add the matrix contribution to diagonal and off diagonal term (matrix is symmetric)
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        //off diagonal
                        res(i1*DIM+d0,i2*DIM+d0)=-sub[d0*DIM+d0];
                        //off diagonal symmetric
                        res(i2*DIM+d0,i1*DIM+d0)=-sub[d0*DIM+d0];
                        //diagonal 1
                        res(i1*DIM+d0,i1*DIM+d0)+=sub[d0*DIM+d0];
                        //diagonal 2
                        res(i2*DIM+d0,i2*DIM+d0)+=sub[d0*DIM+d0];
                        for (unsigned int d1=0;d1<d0;d1++) {
                            // sub[d0*DIM+d1]=sub[d1*DIM+d0];

                            //off diagonal
                            res(i1*DIM+d0,i2*DIM+d1)=-sub[d1*DIM+d0];
                            res(i1*DIM+d1,i2*DIM+d0)=-sub[d1*DIM+d0];

                            //off diagonal symmetric
                            res(i2*DIM+d0,i1*DIM+d1)=-sub[d1*DIM+d0];
                            res(i2*DIM+d1,i1*DIM+d0)=-sub[d1*DIM+d0];

                            //diagonal 1
                            res(i1*DIM+d0,i1*DIM+d1)=sub[d1*DIM+d0];
                            res(i1*DIM+d1,i1*DIM+d0)=sub[d1*DIM+d0];
                            //diagonal 2
                            res(i2*DIM+d0,i2*DIM+d1)=sub[d1*DIM+d0];
                            res(i2*DIM+d1,i2*DIM+d0)=sub[d1*DIM+d0];
                        }
                    }
                }
            }
            res1=res1/2.0;
            return res;
        } else {
            throw std::runtime_error("Not implemented!");
        }

    }
    virtual void init_pbc(const Eigen::Matrix<double,DIM,DIM> &t){
        T=t;
        Tinv=t.inverse();
    }

    bool check_forces(const Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> > & x, const double & dx_over_x, const double & max_error=0.001 ) {
        Eigen::Matrix<double,N*DIM,1> force=deriv(x),x_;
        double res=true;
        for (unsigned int i=0;i<N*DIM;i++) {
            double dvdr=0.0;
            x_=x;
            x_(i)=x_(i)+x_(i)*dx_over_x;
            double vp=operator ()(x_);
            x_(i)=x_(i)-x_(i)*dx_over_x;
            double vm=operator ()(x_);
            dvdr=(vp-vm)/(2*x(i)*dx_over_x);
            std::cerr << "Index "<< i<< ": \tnumerical = " << dvdr << " calculated = "<< force(i)<<"\n";
            if (fabs((dvdr-force(i)))/force(i)>max_error) {
                res=false;
            }
        }
        return res;
    }
    bool check_hessian_forces(const Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> > & x, const double & dx_over_x, const double & max_error=0.001 ) {
        double res=false;
        Eigen::Matrix<double,N*DIM,1> force,x_;
        Eigen::Matrix<double,N*DIM,N*DIM> H=hessian_deriv(x,force);
        for (unsigned int i=0;i<N*DIM;i++) {
            {
                double dvdr=0.0;
                x_=x;
                x_(i)=x_(i)+x_(i)*dx_over_x;
                double vp=operator ()(x_);
                x_(i)=x_(i)-x_(i)*dx_over_x;
                double vm=operator ()(x_);
                dvdr=(vp-vm)/(2*x(i)*dx_over_x);
                std::cerr << "Force "<< i<< ": \t numerical = " << dvdr << " calculated = "<< force(i)<<"\n";
                if (fabs((dvdr-force(i)))/force(i)>max_error) {
                    res=false;
                }
            }
            for (unsigned int j=0;j<N*DIM;j++) {
                double dp,dm,d1p,d1m;
                x_=x;
                x_(j)=x_(j)+x_(j)*dx_over_x;
                x_(i)=x_(i)+x_(i)*dx_over_x;
                dp=operator ()(x_);
                x_=x;
                x_(j)=x_(j)+x_(j)*dx_over_x;
                x_(i)=x_(i)-x_(i)*dx_over_x;
                dm=operator ()(x_);
                d1p=(dp-dm)/(2*x(i)*dx_over_x);

                x_=x;
                x_(j)=x_(j)-x_(j)*dx_over_x;
                x_(i)=x_(i)+x_(i)*dx_over_x;
                dp=operator ()(x_);
                x_=x;
                x_(j)=x_(j)-x_(j)*dx_over_x;
                x_(i)=x_(i)-x_(i)*dx_over_x;
                dm=operator ()(x_);
                d1m=(dp-dm)/(2*x(i)*dx_over_x);

                double d2=(d1p-d1m)/(2*x(j)*dx_over_x);

                std::cerr << i << " "<<j << " " << d2 << " " << H(i,j) << "\n";
                if (fabs((d2-H(i,j)))/H(i,j)>max_error) {
                    res=false;
                }
            }

        }

        std::cout << "\n\n\n"<<H<<"\n";
        return res;
    }

    void pbc_wrap(Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > x ,double L=0.0) {
        if (L==0.0) L=T(0,0);
        x=x-L*Eigen::floor(x.array()/L).matrix();
    }

private:
    constexpr bool want_pbc()       {return (flags & 1) == 1;}
    constexpr bool has_deriv2()     {return (flags & 2) == 2;}
    constexpr bool newton_forces()  {return (flags & 4) == 4;}
    constexpr bool orthorombic_box() {return (flags & 8) == 8;}
    Eigen::Matrix<double,DIM,DIM> T,Tinv;
    template <typename D1,typename D2>
    inline
    Eigen::Matrix<typename D1::Scalar,DIM,1>  pbc(const Eigen::MatrixBase<D1> &x1, const Eigen::MatrixBase<D2> &x0,
                       typename D1::Scalar & r2min) {
        Eigen::Matrix<typename D1::Scalar,DIM,1> u1,u0,dx;
        u1=Tinv*(x1-x0);
        u0=u1-Eigen::round(u1.array()).matrix();
        dx=T*u0;
        r2min=dx.squaredNorm();
        if (! orthorombic_box()){
            for (unsigned int i=0;i<DIM;i++) {
                for (int i0=-1;i0<2;i0+=2){
                    Eigen::Matrix<typename D1::Scalar,DIM,1> dxn=dx+double(i0)*T.col(i);
                    typename D1::Scalar r2=dxn.squaredNorm();
                    if (r2<r2min){
                        r2min=r2;
                        dx=dxn;
                    }
                }
            }
        }
        return dx;
    }
    virtual double pair       (const double & r2)=0;
    virtual double pair_deriv_r2 (const double & r2)=0;
    virtual double pair_deriv2_r2(const double & r2)=0;
};

template <unsigned int N,unsigned int DIM>
class LJPair : public MultiPair<N,DIM,1 | 2 | 8 > {
protected:
    virtual double pair       (const double & r2) override {
        double r6=r2*r2*r2;
        return 1.0/(r6*r6)-1.0/r6;
    }
    virtual double pair_deriv_r2 (const double & r2) override {
        double r6=r2*r2*r2;
        return -6.0/(r6*r6*r2)+3.0/(r6*r2);
    }

    virtual double pair_deriv2_r2(const double & r2) override{
        double r4=r2*r2;
        double r6=r4*r2;
        double r12=r6*r6;
        return 7.0*6.0/(r12*r4)-4.0*3.0/(r6*r4);
    }

};

template <unsigned int N,unsigned int DIM,unsigned int FLAGS>
class Integrator {
public:
    Integrator (MultiPair<N,DIM,FLAGS> * p) : p(p) {

    }

    virtual void step(Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > pos ) {
        std::cerr << "Error: not implemented\n";
        abort();
    }
    virtual void step(Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > pos,Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > vel  ){
        std::cerr << "Error: not implemented\n";
        abort();
    }

    void dump(std::ostream & out, const Eigen::Ref< const Eigen::Matrix<double,N*DIM,1> > & x) {

        out << N <<"\n\n";
        for (unsigned int i=0;i<N;i++) {
            out  <<"1 ";
             for (unsigned int j=0;j<DIM;j++)
                out << x(i*DIM+j) << " ";
             out << "\n";
        }
       // out << "\n;

    }


protected:
    MultiPair<N,DIM,FLAGS> *p;
    double energy;
    Eigen::Matrix<double,N*DIM,1> pos_m;
    Eigen::Matrix<double,N*DIM,1> deriv_m;
    Eigen::Matrix<double,N*DIM,N*DIM> hessian_m;
};
template <unsigned int N, unsigned int DIM,unsigned int FLAGS>
class IntegratorAcceleratedLangevin : public Integrator<N,DIM,FLAGS> {
public:
    static double regularizer(const double & e,double * const & p) {
        /* eq. (45) of Mazzola, Sorella, "Accelerated molecular dynamics for ab-initio Electronic Simulations""
         * e    --> lambda
         * p[0] --> epsilon
         * p[1] --> tau
         * p[2] --> delta
        */
        return 1.0/(1.0+exp((e-p[0])/p[1]))/p[2]+e*(1.0-1.0/(1.0+exp((e-p[0])/p[1])));
    }
    IntegratorAcceleratedLangevin (MultiPair<N,DIM,FLAGS> *p,Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > pos, double delta,double T) : Integrator<N,DIM,FLAGS>(p), delta(delta), T(T) {
        pos_m=pos;
        c=sqrt(2*T*delta);
    }

    void set_T(double T_) {
        T=T_;
        c=sqrt(2*T*delta);
    }

    virtual void step(Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > pos, bool accelerated=true) final {
        Eigen::Matrix<double,N*DIM,1> z,pos_p;
        Eigen::Matrix<double,N*DIM,N*DIM,1> H,H_inv;
        double d=1.0e-1;
        double regularizer_parameters[]={1.0,2/d,d};

        if (accelerated){
            H=p->hessian_deriv(pos,deriv_m);
            //regularize the matrix: diagonalize it
            Eigen::SelfAdjointEigenSolver< Eigen::Matrix<double,N*DIM,N*DIM,1> > eigensolver;
            eigensolver.compute(H);
            auto eigenvalue=eigensolver.eigenvalues();
            auto eigenvectors=eigensolver.eigenvectors();
            //regularize the eigenvalues with regularizer function
            eigenvalue=eigenvalue.unaryExpr(std::bind(&regularizer,std::placeholders::_1,regularizer_parameters));
            H=eigenvectors.transpose()*eigenvalue.asDiagonal()*eigenvectors;
            H_inv=H.inverse();
        } else{
            deriv_m=p->deriv(pos);
        }
        /*
        auto eigenvalues=H.eigenvalues();
        std::cout << eigenvalues<<"\n";
        //if (eigenvalues.col(0)[0]<0){
        //    H=H+eigenvalues(0)*1.1*Eigen::Matrix<double,N*DIM,N*DIM>::Identity();
        //}
        */
        //mette in z un vettore di variabili normali distribuite come N(0,1)
        for (unsigned int i=0;i<N*DIM;i++) {
            z(i)=normal_gauss();
        }

        //trasforma le variabili secondo la matrice di covarianze
        if (accelerated)
            z=H_inv*z;

        if (accelerated){
            pos_p=pos-c*z-H_inv*deriv_m-
                    H_inv*(hessian_m-H)*(pos_m-pos)/2.0;
            hessian_m=H;
        }
        else
            pos_p=pos-c*z -delta*deriv_m;

        pos_m=pos;
        pos=pos_p;
    }

    virtual void step(Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > pos,Eigen::Ref < Eigen::Matrix<double,N*DIM,1> > vel  ) {
        std::cerr << "Warning: called method for second order dynamics, but this is first order\n";
        step(pos);
    }

    void init_global_rnd(unsigned int seed,unsigned int thermalization_steps=10000) {
        set_SHR3_jsr((seed+1)*123);
        for (unsigned int i=0;i<thermalization_steps;i++) {
            rnd_shr3();
        }
        init_cmwc4096();
        for (unsigned int i=0;i<thermalization_steps;i++) {
            cmwc4096();
        }
    }

private:
    double delta,T,c;
    using Integrator<N,DIM,FLAGS>::p;
    using Integrator<N,DIM,FLAGS>::pos_m;
    using Integrator<N,DIM,FLAGS>::deriv_m;
    using Integrator<N,DIM,FLAGS>::hessian_m;

};

int main() {

    bool accelerated=false, test_hessian=false;

    std::ifstream inputjs("input.json");
    json js;
    inputjs >> js;

    /// test con una forma quadratica:

    Eigen::Matrix4d m=Eigen::Matrix4d::Random();
    m = Eigen::Matrix4d::Random();
    Eigen::Matrix<double,4,1> xq=2*Eigen::Matrix<double,4,1>::Random(), b=Eigen::Matrix<double,4,1>::Random();

    QuadraticForm<4> testq;
    testq.set_A_b(0.5*(m.transpose()+m)+Eigen::Matrix4d::Identity()*4,b);
    ParabolaLineMinimization <Eigen::Matrix<double,4,1>,double,QuadraticForm<4> > lineMinq(5,4,20,3);
    Cg<Eigen::Matrix<double,4,1>,double,QuadraticForm<4> > testcgq(testq,xq,testq(xq),lineMinq,8);
    std::cout << 0 << " " << testq.get_solution_distance(testcgq.get_x())<< "\n";
    for (unsigned int i=0;i<4;i++){
        testcgq.iteration();
        std::cout << i+1 << " " << testq.get_solution_distance(testcgq.get_x())<< "\n";
    }

    std::cout << "======================\n";

    LJPair<NATOMS,3> test;
    Eigen::Matrix3d cel;

    if (js.count("cell_size")==0) {
        std::cerr << "Errore: impossibile trovare l'elemento \"cell_size\"\n";
        return -1;
    }

    if (js.count("test_hessian")!=0) {
        test_hessian=js["test_hessian"];
        std::cerr << "test_hessian: " << (test_hessian?"true":"false") << "\n";
    }



    if (js.count("dynamics")==0) {
        std::cerr << "Errore: impossibile trovare l'elemento \"dynamics\"\n";
        return -1;
    } else {
        if (js["dynamics"].count("nsteps")==0) {
            std::cerr << "Errore: impossibile trovare l'elemento \"dynamics\":\"nsteps\"\n";
            return -1;
        }
        if (js["dynamics"].count("output")==0) {
            std::cerr << "Errore: impossibile trovare l'elemento \"dynamics\":\"output\"\n";
            return -1;
        }
        if (js["dynamics"].count("T")==0) {
            std::cerr << "Errore: impossibile trovare l'elemento \"dynamics\":\"T\"\n";
            return -1;
        }
        if (js["dynamics"].count("dt")==0) {
            std::cerr << "Errore: impossibile trovare l'elemento \"dynamics\":\"dt\"\n";
            return -1;
        }
        if (js["dynamics"].count("accelerated")>0) {
            std::cerr << "Accelerated dynamics: ";
            accelerated=js["dynamics"]["accelerated"];
            std::cerr << (accelerated?"true":"false")<<"\n";
        }
    }
    if (js.count("minimization")==0) {
        std::cerr << "Errore: impossibile trovare l'elemento \"minimization\"\n";
        return -1;
        if (js["minimization"].count("nsteps")==0) {
            std::cerr << "Errore: impossibile trovare l'elemento \"minimization\":\"nsteps\"\n";
            return -1;
        }
    }

    unsigned int nsteps_d=js["dynamics"]["nsteps"];
    unsigned int nsteps_cg=js["minimization"]["nsteps"];
    unsigned int dump_mod=1;
    double cell_size=js["cell_size"];
    double temperature=js["dynamics"]["T"];
    double dt=js["dynamics"]["dt"];
    double Tfinal;
    std::string outname=js["dynamics"]["output"];
    if (js["dynamics"].count("print")==1){
        dump_mod=js["dynamics"]["print"];
    }

    if (js["dynamics"].count("Tfinal")==0) {
        Tfinal=temperature;
    } else {
        Tfinal=js["dynamics"]["Tfinal"];
    }

    cel << cell_size , 0.0 , 0.0
        , 0.0 , cell_size , 0.0
        , 0.0 , 0.0 , cell_size;
    test.init_pbc( cel);
    Eigen::Matrix<double,NATOMS*3,1> x=Eigen::Matrix<double,NATOMS*3,1>::Random()*cell_size;

    ParabolaLineMinimization <Eigen::Matrix<double,NATOMS*3,1>,double,LJPair<NATOMS,3> > lineMin(0.02,0.1,4,3);
    Cg<Eigen::Matrix<double,NATOMS*3,1>,double,LJPair<NATOMS,3> > testcg(test,x,test(x),lineMin,8);
    for (unsigned int i=0;i<nsteps_cg;i++) {
        if (i%100==0)
            std::cout << i << " " << testcg.get_fx()<< "\n";
        if (!testcg.iteration())
            break;
    }
    std::cout <<  "Final: " << testcg.get_fx()<< "\n";
    x=testcg.get_x();

    if (test_hessian){
        std::cout << "Test delle forze: " << test.check_forces(x,0.0001,0.001)<<"\n";
        std::cout << "Test dell'hessiana: " << test.check_hessian_forces(x,0.0001,0.001)<<"\n";
    }

    std::cout << "\n\n\nInizio dinamica\n\n\n";
    std::ofstream output(outname,std::ios_base::app);
    IntegratorAcceleratedLangevin<NATOMS,3,11> firstOrderAcceleratedLangevin(&test,x,dt,temperature);
    //forse ok in c++17
//    IntegratorAcceleratedLangevin<> firstOrderAcceleratedLangevin(&test,x,0.001,0.5);
    firstOrderAcceleratedLangevin.init_global_rnd(67857);
//    test.pbc_wrap(x);
    firstOrderAcceleratedLangevin.dump(output,x);
    for (unsigned int istep=0;istep<nsteps_d;istep++) {
        firstOrderAcceleratedLangevin.set_T(temperature+ (Tfinal-temperature)*double(istep)/double(nsteps_d-1));
        firstOrderAcceleratedLangevin.step(x,accelerated);
        if (istep%dump_mod==0)
            firstOrderAcceleratedLangevin.dump(output,x);
        std::cout <<istep<<" " << test(x) << "\n";
    }
    std::cout << "Finito!\nVMD:\n\n";
    std::cout <<"mol delete 0\n\
mol addrep 0\n\
display resetview\n\
mol new {/home/bertossa/analisi/build-analisi-Desktop-Minimum Size Release/libcg/out.xyz} type {xyz} first 0 last -1 step 1 waitfor -1\n"
    << "set cell [pbc set { "<< cell_size <<" " <<cell_size<<" " <<cell_size<< " } -all]\npbc wrap -all\npbc box\nmol modstyle 0 0 VDW 0.100000 12.000000\n";

}
