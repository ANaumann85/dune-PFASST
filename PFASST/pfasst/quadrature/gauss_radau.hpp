/**
 * @file pfasst/quadrature/gauss_radau.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__QUADRATURE__GAUSS_RADAU_HPP_
#define _PFASST__QUADRATURE__GAUSS_RADAU_HPP_

#include "pfasst/quadrature/quadrature.hpp"


namespace pfasst
{
  namespace quadrature
  {
    /**
     * Quadrature handler for Gauss-Radau quadrature.
     *
     * @tparam scalar precision of quadrature (i.e. `double`)
     *
     * @ingroup Quadrature
     * @since v0.3.0
     */
    template<typename precision = pfasst::time_precision>
    class GaussRadau
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = false;
        static const bool RIGHT_IS_NODE = true;
        //! @}

      public:
        //! @{
        /**
         * @throws invalid_argument if less than two nodes are requested
         */
        explicit GaussRadau(const size_t num_nodes);
        GaussRadau() = default;
        virtual ~GaussRadau() = default;
        //! @}

        //! @{
        virtual bool left_is_node() const override;
        virtual bool right_is_node() const override;
        //! @}
        
        //! @{
        virtual string print_summary() const override;
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override;
        virtual void compute_weights() override;
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#include "pfasst/quadrature/gauss_radau_impl.hpp"

#endif  // _PFASST__QUADRATURE__GAUSS_RADAU_HPP_
