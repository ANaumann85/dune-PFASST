#ifndef _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_

//#include <dune/istl/matrixmarket.hh>
#include <dune/common/densematrix.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>

#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/pq1nodalbasis.hh>
#include <memory>
#include <type_traits>

using std::shared_ptr;
using std::vector;

#include <vector>

#include <pfasst/sweeper/FE_imex.hpp>
#include <pfasst/contrib/fft.hpp>
#include "fe_manager.hpp"

//using namespace Dune;

namespace pfasst
{
  namespace examples
  {
    namespace heat_FE
    {
        template<
                class EncapsulationTraits,
                size_t base_order,
                size_t dim,
                class... Ts
        >
        struct dune_sweeper_traits : public sweeper_traits<EncapsulationTraits,  Ts...>
        {
            //! type of the Encapsulation traits
            using encap_traits = EncapsulationTraits;
            //! type of the user data encapsulation
            using encap_t = encap::Encapsulation<EncapsulationTraits>;
            //! type of the temporal domain
            using time_t = typename encap_traits::time_t;
            //! type of the spacial domain
            using spatial_t = typename encap_traits::spatial_t;

            //static constexpr size_t BASE_ORDER = base_order;
            //static constexpr size_t DIM = dim;
        };

      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class Heat2d_FE
        : public pfasst::IMEX<SweeperTrait, Enabled>
      {
        static_assert(std::is_same<
                        typename SweeperTrait::encap_t::traits::dim_t,
                        std::integral_constant<size_t, 1>
                      >::value,
                      "Heat2d_FE Sweeper requires 2D data structures");

        public:
          using traits = SweeperTrait;

          static void init_opts();

          int finer;
	  static const unsigned                          dim = 2;  
	  
	  typedef Dune::YaspGrid<dim> GridType;
	  struct MovingFunction : Dune::AnalyticalCoordFunction< double, dim, dim,  MovingFunction>
	  {
	    typedef Dune::AnalyticalCoordFunction< double, dim, dim,  MovingFunction> super;
	    typedef typename super::DomainVector  DomainVector;
	    typedef typename super::RangeVector   RangeVector;
	    
	    std::array<double, dim > dx;
	    
	    MovingFunction()
	    {
	      for(unsigned i(0); i < dim; ++i)
		dx[i] = 0.0;
	    }
	    
	    void evaluate ( const DomainVector &x, RangeVector &y ) const
	    {
	      for(unsigned i(0); i < dim ; ++i)
		y[i] = x[i]+dx[i];
	    }
	  };
        
	  typedef Dune::YaspGrid<dim,Dune::EquidistantOffsetCoordinates<double,dim> > HostGridType_MV;
	  typedef Dune::GeometryGrid<HostGridType_MV, MovingFunction >                GridType_MV;
	  typedef typename GridType_MV::LeafGridView                                  GridView_MV;
      private:
          using spatial_t = typename traits::spatial_t;

          typename traits::time_t                        _t0{0.0};
          spatial_t                                      _nu{1.0e-3};
          pfasst::contrib::FFT<typename traits::encap_t> _fft;
          vector<vector<spatial_t>>                      _lap;
	 


          typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;
          typedef Dune::BlockVector<Dune::FieldVector<double,1> > VectorType;

          //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dune;
          //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > A_dune;
	  //VectorType b_dune; //the constant source for testing purposes
	  
          //Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > M_dtA_dune;

      

          std::shared_ptr<GridType> grid;
          std::shared_ptr<GridType> grid_coarse;
	  



          typedef GridType::LeafGridView GridView;


          using BasisFunction = Dune::Functions::PQ1NodalBasis<GridView>;

          std::shared_ptr<BasisFunction> basis;
          std::shared_ptr<BasisFunction> basis_coase;
	  
	  //moving information
	  const Dune::FieldVector< double, dim > lower_mv, upper_mv;
	  mutable MovingFunction mf;	  
	  HostGridType_MV hgtmv;
	  std::shared_ptr< GridType_MV> grid_mv;
	  double v0, alpha, sourceVal, vel;
	  void fastGrid(double t, const VectorType& yIn, VectorType& out) const;
	  void fastRect(double t, const VectorType& yIn, VectorType& out) const;
	  template<typename Flux >
	  void fastBoundary(const VectorType& yIn, const Flux& flux, VectorType& out) const;


        protected:
#if 1
          virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_expl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;
#endif

          virtual shared_ptr<typename SweeperTrait::encap_t>
          evaluate_rhs_impl(const typename SweeperTrait::time_t& t,
                            const shared_ptr<typename SweeperTrait::encap_t> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_t> f,
                                      shared_ptr<typename SweeperTrait::encap_t> u,
                                      const typename SweeperTrait::time_t& t,
                                      const typename SweeperTrait::time_t& dt,
                                      const shared_ptr<typename SweeperTrait::encap_t> rhs) override;

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_error(const typename SweeperTrait::time_t& t);

          virtual vector<shared_ptr<typename SweeperTrait::encap_t>>
          compute_relative_error(const vector<shared_ptr<typename SweeperTrait::encap_t>>& error,
                                 const typename SweeperTrait::time_t& t);

        public:
          explicit Heat2d_FE(const size_t nelements, const size_t basisorder, const int finer);


//           Heat2d_FE(const Heat_FE<SweeperTrait, Enabled>& other)= default;
//           Heat2d_FE(Heat_FE<SweeperTrait, Enabled>&& other)= default;
//           virtual ~Heat2d_FE() = default;
//           Heat_FE<SweeperTrait, Enabled>& operator=(const Heat_FE<SweeperTrait, Enabled>& other) = default;
//           Heat_FE<SweeperTrait, Enabled>& operator=(Heat_FE<SweeperTrait, Enabled>&& other) = default;

          virtual void set_options() override;

          template<typename Basis>
          void assemble(Basis &basis);

          void assemble();

          virtual shared_ptr<typename SweeperTrait::encap_t> exact(const typename SweeperTrait::time_t& t);
	  shared_ptr<typename SweeperTrait::encap_t> initial();
          virtual void post_step() override;

          virtual bool converged(const bool pre_check) override;
          virtual bool converged() override;

          size_t get_num_dofs() const;

          shared_ptr<GridType> get_grid() const;
	  static shared_ptr<GridType> createGrid(unsigned nElements);
      };
    }  // ::pfasst::examples::heat_FE
  }  // ::pfasst::examples
}  // ::pfasst

#include "FE_heat2d_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD2D__HEAD2D_SWEEPER_HPP_
