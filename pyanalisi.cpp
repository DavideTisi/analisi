#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "traiettoria.h"
#include "readlog.h"
#include "correlatorespaziale.h"
#include "gofrt.h"
#include "greenkuboNcomponentionicfluid.h"

namespace py = pybind11;

PYBIND11_MODULE(pyanalisi,m) {
    py::class_<Traiettoria>(m,"Traj")
            .def(py::init<std::string>(),R"begend(
                 Parameters
                 ----------
                 name pf the binary file to open
)begend")
            .def("setWrapPbc",&Traiettoria::set_pbc_wrap,R"begend(
                 wrap all the atomic coordinate inside the simulation box, red from the binary file

                 Parameters
                 ----------
                 True/False
)begend")
            .def("getNtimesteps",&Traiettoria::get_ntimesteps,R"begend(
                 returns estimated number of timesteps from the file size
)begend")
            .def("setAccessWindowSize",[](Traiettoria & t,int ts) {return (int) t.imposta_dimensione_finestra_accesso(ts);},R"begend(
                 sets the size of the read block. Must fit in memory
                 
                 Parameters
                 ----------
                 int -> size of the block
                
                 Returns
                 ----------
                 1 -> success
                 0/-1 -> failure
)begend")
            .def("setAccessStart",[](Traiettoria & t,int ts) {return (int) t.imposta_inizio_accesso(ts);},R"begend(
                 sets the first timestep to read, and read the full block

                 Parameters
                 ----------
                 int -> first timestep to read
)begend");

    using RL = ReadLog<double>;
    py::class_<RL>(m,"ReadLog")
            .def(py::init<std::string,Traiettoria*,unsigned int, unsigned int, unsigned int, std::vector<std::string> >() ,R"begend(
                 perform an advanced read of a column formatted textfile with headers. Support multithreaded reading of chunks of lines.

                 Parameters
                 ----------
                 string -> filename
                 Traj -> Traj instance, to eventually calculate some currents
                 unsigned int ->
                 unsigned int ->
                 unsigned int ->
                 string list -> 
)begend")
            .def("getNtimesteps", &RL::n_timestep,R"begend(
                 return number of timesteps read
)begend");

    py::class_<CorrelatoreSpaziale>(m,"CurrentCalc", py::buffer_protocol())
            .def(py::init<Traiettoria*,std::vector< std::array<double,3> >,double,unsigned int,unsigned int,bool>(),R"begend(
                 calculates  \dot \tilde e(k) / |k|  (can be useful to define some new current)
)begend")
            .def("reset",&CorrelatoreSpaziale::reset,R"begend()begend")
            .def("calculate",&CorrelatoreSpaziale::calcola,R"begend()begend")
            .def("printRes",[](CorrelatoreSpaziale & c){
                c.print(std::cout);
            },R"begend()begend")
            .def("__enter__",[](CorrelatoreSpaziale & c) -> CorrelatoreSpaziale & { return c;} )
            .def("__exit__",[](CorrelatoreSpaziale & c, py::object * exc_type, py::object * exc_value, py::object * traceback){})
            .def_buffer([](CorrelatoreSpaziale & c) -> py::buffer_info{
                return py::buffer_info(
                c.accesso_lista(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                c.get_shape().size(),                                      /* Number of dimensions */
                c.get_shape(),                 /* Buffer dimensions */
                c.get_stride()
        );
            });

    py::class_<Gofrt<double> >(m,"Gofrt", py::buffer_protocol())
            .def(py::init<Traiettoria*,double,double,unsigned int, unsigned int, unsigned int,unsigned int, bool>(),R"begend(
                 calculates g(r) and, in general, g(r,t). Interesting dynamical property
)begend")
            .def("reset",& Gofrt<double>::reset,R"begend()begend")
            .def("getNumberOfExtraTimestepsNeeded",&Gofrt<double>::numeroTimestepsOltreFineBlocco ,R"begend()begend")
            .def("calculate", & Gofrt<double>::calcola,R"begend()begend")
            .def_buffer([](Gofrt<double> & g) -> py::buffer_info {
                return py::buffer_info(
                g.accesso_lista(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                g.get_shape().size(),                                      /* Number of dimensions */
                g.get_shape(),                 /* Buffer dimensions */
                g.get_stride());
            });

    using GK = GreenKuboNComponentIonicFluid<double,double>;
    py::class_<GK>(m,"GreenKubo", py::buffer_protocol())
            .def(py::init<RL*, std::string, unsigned int, std::vector<std::string>, bool, unsigned int, unsigned int, bool, unsigned int, unsigned int, bool, unsigned int, unsigned int>(),R"begend()begend")
            .def("getNumberOfExtraTimestepsNeeded",&GK::numeroTimestepsOltreFineBlocco,R"begend(
                perform an efficient, multithreaded Green-Kubo calculations of the transport coefficient of a condensed matter system, given the currents.
)begend")
            .def("reset",&GK::reset,R"begend()begend")
            .def("calculate",&GK::calcola,R"begend()begend")
            .def_buffer([](GK & g) -> py::buffer_info {
                return py::buffer_info(
                g.accesso_lista(),                               /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                g.get_shape().size(),                                      /* Number of dimensions */
                g.get_shape(),                 /* Buffer dimensions */
                g.get_stride());
            });


}
