#ifndef CENTERDIFF_H
#define CENTERDIFF_H

#include "calcolamultithread.h"
#include "operazionisulista.h"

class CenterDiff : public CalcolaMultiThread, public OperazioniSuLista<CenterDiff,double>
{
public:
    CenterDiff(Traiettoria *t, unsigned int nthreads=0, unsigned int skip=1, unsigned int nit=1,bool sum_first_two_and_ignore_vz=false,bool sum_1=false);
    virtual unsigned int numeroTimestepsOltreFineBlocco(unsigned int n_b) final {return 0;}
    virtual void reset(const unsigned int numeroTimestepsPerBlocco) final;
    virtual void calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int & ith) noexcept final;
    virtual std::vector<ssize_t> get_shape() const final {return {lunghezza_lista/3/nit/3,nit,3,3}; }
    virtual std::vector<ssize_t> get_stride() const final {return { static_cast<long> (sizeof (double)*3*3*nit),sizeof (double)*3*3,sizeof (double)*3,sizeof (double)};}
    void set_starting_center(const std::array<double,3*3> & s) {starting_center=s;}
    ~CenterDiff();
private:
    Traiettoria *t;
    using OperazioniSuLista<CenterDiff,double>::lista;
    using OperazioniSuLista<CenterDiff,double>::lunghezza_lista;
    unsigned int nit,lista_alloc;
    bool sum_first_two_and_ignore_vz,sum_1;
    std::array<double,3*3> starting_center;

};

#endif // CENTERDIFF_H
