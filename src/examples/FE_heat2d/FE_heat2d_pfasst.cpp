//#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>


#include <dune/istl/matrixmarket.hh>
#include <memory>
#include <stdexcept>
using std::shared_ptr;

#include <mpi.h>
#include "dune_includes"
#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/controller/two_level_pfasst.hpp>

#include <pfasst/comm/mpi_p2p.hpp>
//#include <pfasst/contrib/spectral_transfer.hpp>
#include "FE_heat2d_sweeper.hpp"
#include "../../datatypes/dune_vec.hpp"
#include "spectral_transfer.hpp"

#include <dune/common/timer.hh>

//////////////////////////////////////////////////////////////////////////////////////
//
// Compiletimeparameter
//
//////////////////////////////////////////////////////////////////////////////////////

const size_t DIM = 2;            //Räumliche Dimension des Rechengebiets

const size_t BASIS_ORDER = BASE_ORDER;    //maximale Ordnung der Lagrange Basisfunktionen

//////////////////////////////////////////////////////////////////////////////////////


using encap_traits_t = pfasst::encap::dune_vec_encap_traits<double, double, 1>;
using pfasst::encap::DuneEncapsulation;


using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::contrib::SpectralTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommunicatorType;

using pfasst::examples::heat_FE::Heat2d_FE;

typedef DuneEncapsulation<double, double, 1>                     EncapType;


typedef Heat2d_FE<pfasst::examples::heat_FE::dune_sweeper_traits<encap_traits_t, BASIS_ORDER, DIM>> SweeperType;


typedef pfasst::transfer_traits<SweeperType, SweeperType, 1>       TransferTraits;
typedef SpectralTransfer<TransferTraits>                           TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
      void run_pfasst(const size_t nelements, const size_t basisorder, const size_t dim, const size_t& nnodes, const pfasst::quadrature::QuadratureType& quad_type,
                      const double& t_0, unsigned nSteps, const double& t_end, const size_t& niter)
      {
        using pfasst::config::get_value;
        unsigned nCoarse = get_value<unsigned>("--nCoarse",1);
        unsigned nFine = get_value<unsigned>("--nFine",1);
        typedef SweeperType::GridType GridType;
        int my_rank, num_pro;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank );
        MPI_Comm_size(MPI_COMM_WORLD, &num_pro );

        TwoLevelPfasst<TransferType, CommunicatorType> pfasst;
        pfasst.communicator() = std::make_shared<CommunicatorType>(MPI_COMM_WORLD);
//         pfasst.grid_builder(nelements);
	
        std::shared_ptr<GridType > mlsdcGrid(SweeperType::createGrid(nelements));
        
        auto FinEl = make_shared<fe_manager>(mlsdcGrid, nelements, 2); 
        //pfasst.grid_builder(mlsdcGrid);

        auto coarse = std::make_shared<SweeperType>(nelements, basisorder, 0);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = std::make_shared<SweeperType>(nelements, basisorder, 1);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        coarse->is_coarse = true;
        fine->is_coarse = false;
        auto transfer = std::make_shared<TransferType>();



        pfasst.add_transfer(transfer);
        pfasst.add_sweeper(coarse, true);
        pfasst.add_sweeper(fine);
        pfasst.set_options();


        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = (t_end-t_0)/nSteps;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        //coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());	
        //fine->initial_state() = fine->exact(pfasst.get_status()->get_time());
	
        coarse->initial_state() = coarse->initial();
        fine->initial_state() = fine->initial();
        Dune::Timer timer;
        
        pfasst.run(nCoarse, nFine);
        std::cout << "solve-pfasst: " << timer.elapsed() << std::endl;
        pfasst.post_run();



        MPI_Barrier(MPI_COMM_WORLD);


        if(my_rank==num_pro-1) {
          std::stringstream fname;
          fname << "heat2dMoving_result_pfasst_nCoarse_" << nCoarse << "_nFine_" << nFine << "_nIter_" << niter << "_nSteps_" << nSteps << "_nCore_" << num_pro;
          auto anfang    = fine->exact(0)->data();
          auto naeherung = fine->get_end_state()->data();
          auto exact     = fine->exact(t_end)->data();
          //         for (int i=0; i< fine->get_end_state()->data().size(); i++){
          //           std::cout << anfang[i] << " " << naeherung[i] << "   " << exact[i] << " "  <<  std::endl;
          //         }
          
          std::cout << "******************************************* " << std::endl;
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
            vtkWriter.addVertexData(z, "initial_data");
            
            vtkWriter.write("heat2dMoving_result");
          }  else {
            fname << ".csv";
            auto sweeper = fine;
            VectorType x = sweeper->get_end_state()->data();
            std::ofstream file(fname.str().c_str(), std::ios_base::trunc | std::ios_base::out);
            file << std::setprecision(15);
            for(auto& d : x)
              file << d << std::endl;
            file.close();
          }
          std::cout << "******************************************* " << std::endl;
          std::cout << " " << std::endl;
          std::cout << " " << std::endl;
          
          //         fine->get_end_state()->scaled_add(-1.0, fine->exact(t_end));
          //         std::cout << "Fehler: "  << fine->get_end_state()->norm0() << " " << std::endl;
          
          
        }

      }
    }  // ::pfasst::examples::heat_FE
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char* argv[])
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  MPI_Init(&argc, &argv);

  pfasst::init(argc, argv, SweeperType::init_opts);
  pfasst::Status<double>::create_mpi_datatype();
  pfasst::config::read_commandline(argc, argv);


  const size_t nelements = get_value<size_t>("--num_elements", 4); //Anzahl der Elemente pro Dimension
  const size_t nnodes = get_value<size_t>("--num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  
  double t_end = get_value<double>("tend", -1);
  size_t nsteps = get_value<size_t>("--num_steps", 0);
  
  double dt = get_value<double>("dt", 0.1);
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


  pfasst::examples::heat_FE::run_pfasst(nelements, BASIS_ORDER, DIM, nnodes, quad_type, t_0, nsteps, t_end, niter);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
