#include <dune/istl/matrixmarket.hh>
#include <fenv.h>


#include <memory>
#include <stdexcept>
using std::shared_ptr;
#include "dune_includes"

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/two_level_mlsdc.hpp>

#include "FE_heat2d_sweeper.hpp"
//#include <pfasst/encap/dune_vec.hpp>
#include "../../datatypes/dune_vec.hpp"
#include "spectral_transfer.hpp"
//#include <pfasst/contrib/spectral_transfer.hpp>




// #include <dune/grid/io/file/vtk/vtkwriter.hh>
// #include <dune/grid/yaspgrid.hh>
// #include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
// 
// #include <dune/common/densematrix.hh>
// 
// #include <dune/istl/bvector.hh>
// #include <dune/istl/bcrsmatrix.hh>
// #include <dune/istl/multitypeblockmatrix.hh>
// 
// #include <dune/grid/yaspgrid.hh>
// 
// #include <dune/functions/functionspacebases/pqknodalbasis.hh>
// #include <dune/functions/functionspacebases/pq1nodalbasis.hh>
// #include <dune/typetree/utility.hh>
// 
// #include <dune/fufem/assemblers/transferoperatorassembler.hh>
// 
// 
// 
// 
// #include <vector>
// #include <dune/common/function.hh>
// #include <dune/common/bitsetvector.hh>
// //#include <dune/common/indices.hh>
// #include <dune/geometry/quadraturerules.hh>
// 
// #include <dune/grid/yaspgrid.hh>
// #include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
// 
// #include <dune/istl/matrix.hh>
// #include <dune/istl/bcrsmatrix.hh>
// #include <dune/istl/multitypeblockmatrix.hh>
// 
// #include <dune/istl/multitypeblockvector.hh>
// #include <dune/istl/matrixindexset.hh>
// #include <dune/istl/solvers.hh>
// #include <dune/istl/preconditioners.hh>
// 
// 
// #include <dune/functions/functionspacebases/interpolate.hh>
// 
// #include <dune/functions/functionspacebases/taylorhoodbasis.hh>
// #include <dune/functions/functionspacebases/hierarchicvectorwrapper.hh>
// #include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
// #include <dune/functions/gridfunctions/gridviewfunction.hh>
// #include <dune/common/timer.hh>


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;



//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = DIMENSION;            //Räumliche Dimension des Rechengebiets

const size_t BASIS_ORDER = BASE_ORDER;    //maximale Ordnung der Lagrange Basisfunktionen

//////////////////////////////////////////////////////////////////////////////////////





namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      using pfasst::transfer_traits;
      using pfasst::contrib::SpectralTransfer;
      using pfasst::TwoLevelMLSDC;
      using pfasst::quadrature::QuadratureType;

      using sweeper_t = Heat2d_FE<pfasst::sweeper_traits<encap_traits_t>>;
      using transfer_traits_t = pfasst::transfer_traits<sweeper_t, sweeper_t, 1>;
      using transfer_t = SpectralTransfer<transfer_traits_t>;
      using heat_FE_mlsdc_t = TwoLevelMLSDC<transfer_t>;


      void run_mlsdc(const size_t nelements, const size_t basisorder, const size_t DIM, const size_t coarse_factor,
                                           const size_t nnodes, const QuadratureType& quad_type,
                                           const double& t_0, unsigned nSteps, const double& t_end,
                                           const size_t niter) 
      {
        using pfasst::config::get_value;
        typedef sweeper_t::GridType GridType;
        
        unsigned nCoarse = get_value<std::size_t>("--nCoarse",1);
        unsigned nFine = get_value<std::size_t>("--nFine",1);
        auto mlsdc = std::make_shared<heat_FE_mlsdc_t>();

        std::shared_ptr<GridType > mlsdcGrid(sweeper_t::createGrid(nelements));
        auto FinEl = make_shared<fe_manager>(mlsdcGrid, nelements, 2); 

        //mlsdc->grid_builder(mlsdcGrid);


        using pfasst::quadrature::quadrature_factory;

        auto coarse = std::make_shared<sweeper_t>(nelements, basisorder, 0);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<sweeper_t>(nelements, basisorder, 0);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        coarse->is_coarse = true;
        fine->is_coarse = false;

        auto transfer = std::make_shared<transfer_t>();
        transfer->create(FinEl);




        mlsdc->add_sweeper(coarse, true);
        mlsdc->add_sweeper(fine /*, false*/);

        mlsdc->add_transfer(transfer);
        mlsdc->set_options();


        mlsdc->status()->time() = t_0;
        mlsdc->status()->dt() = (t_end-t_0)/nSteps;
        mlsdc->status()->t_end() = t_end;
        mlsdc->status()->max_iterations() = niter;


        mlsdc->setup();


        coarse->initial_state() = coarse->initial(); //coarse->exact(mlsdc->get_status()->get_time());
        fine->initial_state() = fine->initial(); //fine->exact(mlsdc->get_status()->get_time());

//         for (int i=0; i< fine->initial_state()->data().size(); i++){
//           std::cout << "Anfangswerte feiner Sweeper: " << " " << fine->initial_state()->data()[i] << std::endl;
//         }

//         std::cout  <<  std::endl;

//         for (int i=0; i< coarse->initial_state()->data().size(); i++){
//           std::cout << "Anfangswerte grober Sweeper: " << " " << coarse->initial_state()->data()[i] <<  std::endl;
//         }




        Dune::Timer timer;
        mlsdc->run(nCoarse, nFine);
        std::cout << "solve-mlscd: " << timer.elapsed() << std::endl;


        mlsdc->post_run();


        std::cout <<  "fein" << std::endl;
        auto naeherung = fine->get_end_state()->data();
        auto exact     = fine->exact(t_end)->data();
#if 0
        for (int i=0; i< fine->get_end_state()->data().size(); i++){
          std::cout << fine->exact(0)->data()[i] << " " << naeherung[i] << "   " << exact[i] << std::endl;
        }


        std::cout << "******************************************* " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << " " <<  std::endl ;
        std::cout << "Fehler: " <<  std::endl ;
        fine->states()[fine->get_states().size()-1]->scaled_add(-1.0 , fine->exact(t_end));
        std::cout << fine->states()[fine->get_states().size()-1]->norm0()<<  std::endl ;
        std::cout << "******************************************* " <<  std::endl ;
#endif


        /*ofstream f;
        stringstream ss;
        ss << nelements;
        string s = "neu_mlsdc/" + ss.str() + ".dat";
        f.open(s, ios::app | std::ios::out );
        f << nelements << " " << dt << " "<< fine->states()[fine->get_states().size()-1]->norm0() << endl;
        f.close();*/
        std::stringstream fname;
        fname << "heat2dMoving_result_mlsdc_nCoarse_" << nCoarse << "_nFine_" << nFine << "_nIter_" << niter << "_nSteps_" << nSteps;
        if(BASIS_ORDER==1) {
          auto sweeper = fine;
          auto grid = (*sweeper).get_grid();
          typedef GridType::LeafGridView GridView;
          GridType::LeafGridView gridView = grid->leafGridView();
          Dune::VTKWriter<GridView> vtkWriter(gridView);
          typedef Dune::BlockVector<Dune::FieldVector<double, 1> > VectorType;
          VectorType x = sweeper->get_end_state()->data();
          //VectorType y = sweeper->exact(t_end)->data();
          VectorType z = sweeper->initial_state()->data();
          vtkWriter.addVertexData(x, "fe_solution");
          //vtkWriter.addVertexData(y, "exact_solution");
          //vtkWriter.addVertexData(z, "initial_data");
          
          vtkWriter.write(fname.str());
        } else {
          fname << ".csv";
          auto sweeper = fine;
          VectorType x = sweeper->get_end_state()->data();
          std::ofstream file(fname.str().c_str(), std::ios_base::trunc | std::ios_base::out);
          file << std::setprecision(15);
          for(auto& d : x)
            file << d << std::endl;
          file.close();
        }
        if(BASIS_ORDER==1 && false) {
          auto sweeper = coarse;
          auto grid = (*sweeper).get_grid();
          typedef GridType::LeafGridView GridView;
          GridType::LeafGridView gridView = grid->leafGridView();
          Dune::VTKWriter<GridView> vtkWriter(gridView);
          typedef Dune::BlockVector<Dune::FieldVector<double, 1> > VectorType;
          VectorType x = sweeper->get_end_state()->data();
          //VectorType y = sweeper->exact(t_end)->data();
          VectorType z = sweeper->initial_state()->data();
          vtkWriter.addVertexData(x, "fe_solution");
          //vtkWriter.addVertexData(y, "exact_solution");
          vtkWriter.addVertexData(z, "initial_data");
          
          vtkWriter.write("heat2dMoving_result_mlsdc_coarse");
        }

      }

    }  // ::pfasst::examples::heat_FE
  } // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{

  feenableexcept(FE_INVALID | FE_OVERFLOW);

  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;
  using sweeper_t = pfasst::examples::heat_FE::Heat2d_FE<pfasst::sweeper_traits<encap_traits_t>>;

  pfasst::init(argc, argv, sweeper_t::init_opts);
  pfasst::config::read_commandline(argc, argv);

  const size_t nelements = get_value<size_t>("--num_elements", 2);
  const size_t nnodes = get_value<size_t>("--num_nodes", 3);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 1);
  const QuadratureType quad_type = QuadratureType::Uniform_Right;
  const double t_0 = 0.0;
  double dt = get_value<double>("dt", 0.001);
  double t_end = get_value<double>("--tend", 0.001);
  size_t nsteps = get_value<size_t>("--num_steps", 0);
  if (t_end == -1 && nsteps == 0) {
    ML_CLOG(ERROR, "USER", "Either t_end or num_steps must be specified.");
    throw std::runtime_error("either t_end or num_steps must be specified");
  } else if (t_end != -1 && nsteps != 0) {
    dt = t_end / nsteps;
#if 0
    if (!pfasst::almost_equal(t_0 + nsteps * dt, t_end)) {
      ML_CLOG(ERROR, "USER", "t_0 + nsteps * dt != t_end ("
                          << t_0 << " + " << nsteps << " * " << dt << " = " << (t_0 + nsteps * dt)
                          << " != " << t_end << ")");
      throw std::runtime_error("t_0 + nsteps * dt != t_end");
    }
#endif
  } else if (nsteps != 0) {
    t_end = t_0 + dt * nsteps;
  }
  const size_t niter = get_value<size_t>("--num_iters", 50);

  pfasst::examples::heat_FE::run_mlsdc(nelements, BASIS_ORDER, DIM, coarse_factor, nnodes, quad_type, t_0, nsteps, t_end, niter);
}
#endif 
