#include "dune_vec.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <type_traits>
#include <vector>
using std::shared_ptr;
using std::vector;

#include "pfasst/logging.hpp"
#include "pfasst/util.hpp"


namespace pfasst
{
  namespace encap
  {
    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const typename traits::GFS_& gfs)
      : _data(gfs) //, gfs(gfs)
    {
      this->zero();
    }

    template<class EncapsulationTrait>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::Encapsulation(const typename EncapsulationTrait::data_t& data)
      : Encapsulation<EncapsulationTrait>(data.N())
    {
      this->data() = data;
    }

    template<class EncapsulationTrait>
    Encapsulation<EncapsulationTrait>&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::operator=(const typename EncapsulationTrait::data_t& data)
    {
      this->data() = data;
      return *this;
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::data()
    {
      return this->_data;
    }

    template<class EncapsulationTrait>
    const typename EncapsulationTrait::data_t&
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_data() const
    {
      return this->_data;
    }

    template<class EncapsulationTrait>
    size_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_total_num_dofs() const
    {
      return this->get_data().N();
    }

    template<class EncapsulationTrait>
    std::array<int, EncapsulationTrait::DIM>
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::get_dimwise_num_dofs() const
    {
      std::array<int, EncapsulationTrait::DIM> dimwise_ndofs;
      switch (EncapsulationTrait::DIM) {
        case 1:
          dimwise_ndofs.fill((int)this->get_total_num_dofs()); //setzt 1,2 v 3 dimensionen und quadratischen raum vorraus
          break;
        case 2:
          dimwise_ndofs.fill((int)sqrt(this->get_total_num_dofs()));
          break;
        case 3:
          dimwise_ndofs.fill((int)cbrt(this->get_total_num_dofs()));
          break;
        default:
          ML_CLOG(FATAL, "ENCAP", "unsupported spatial dimension: " << EncapsulationTrait::DIM);
          throw std::runtime_error("unsupported spatial dimension");
      }

      return dimwise_ndofs;
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::zero()
    {
      std::fill(this->data().begin(), this->data().end(), typename EncapsulationTrait::spatial_t(0.0));
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::scaled_add(const typename EncapsulationTrait::time_t& a,
                                   const shared_ptr<Encapsulation<EncapsulationTrait>> y)
    {

      assert(this->get_data().N() == y->data().N());

      std::transform(this->get_data().begin(), this->get_data().end(), y->data().begin(),
                this->data().begin(),
                [a](const typename EncapsulationTrait::spatial_t& xi,
                    const typename EncapsulationTrait::spatial_t& yi) { return xi + a * yi; });

    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spatial_t
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::norm0() const
    {
      return std::abs(*(std::max_element(this->get_data().begin(), this->get_data().end(),
                               [](const typename EncapsulationTrait::spatial_t& a,
                                  const typename EncapsulationTrait::spatial_t& b)
                                 { return std::abs(a) < std::abs(b); })));
    }

    template<class EncapsulationTrait>
    template<class CommT>
    bool
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::probe(shared_ptr<CommT> comm, const int src_rank, const int tag)
    {
      return comm->probe(src_rank, tag);
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::send(shared_ptr<CommT> comm, const int dest_rank,
                              const int tag, const bool blocking)
    {
      ML_CVLOG(2, "ENCAP", "sending data: " << this->get_data());
      if (blocking) {
        //comm->send(this->get_data().data(), this->get_data().size(), dest_rank, tag);
        comm->send(&(this->data()[0][0]), this->size, dest_rank, tag);
      } else {
        //comm->isend(this->get_data().data(), this->get_data().size(), dest_rank, tag);
        comm->isend(&(this->data()[0][0]), this->size, dest_rank, tag);
      }
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::recv(shared_ptr<CommT> comm, const int src_rank,
                              const int tag, const bool blocking)
    {
      if (blocking) {
        //comm->recv(this->data().data(), this->get_data().size(), src_rank, tag);
        comm->recv(&(this->data()[0][0]), this->size, src_rank, tag);
      } else {
        //comm->irecv(this->data().data(),    this->get_data().size(), src_rank, tag);
        comm->irecv(&(this->data()[0][0]), this->size, src_rank, tag);
      }
      ML_CVLOG(2, "ENCAP", "received data: " << this->get_data());
    }

    template<class EncapsulationTrait>
    template<class CommT>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::bcast(shared_ptr<CommT> comm, const int root_rank)
    {
      //comm->bcast(this->data().data(), this->get_data().size(), root_rank);
      comm->bcast(&(this->data()[0][0]), this->size, root_rank);
    }

    template<class EncapsulationTrait>
    void
    Encapsulation<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::log(el::base::type::ostream_t& os) const
    {
      os << "pdeVector ";// << pfasst::join(this->get_data(), ", ");
    }


    template<class EncapsulationTrait>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::EncapsulationFactory(const size_t size)
      : _size(size)
    {}

    template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if< 
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::create(EncapsulationTrait::GFS::Traits::gridView gv, EncapsulationTrait::GFS::Traits::FiniteElementMap fem) const
    {
      return std::make_shared<Encapsulation<EncapsulationTrait>>(gv, fem);
    }

    template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_size(const size_t& size)
    {
      this->_size = size;
    }


    template<class EncapsulationTrait>
    void
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::set_gfs(typename EncapsulationTrait::GFS_& gfs_)
    {
      this->gfs =  std::make_shared<typename EncapsulationTrait::GFS_>(gfs_);
    }



    template<class EncapsulationTrait>
    size_t
    EncapsulationFactory<
      EncapsulationTrait,
      typename std::enable_if<
                 std::is_same<dune_encap_tag, typename EncapsulationTrait::tag_t>::value
               >::type>::size() const
    {
      return this->_size;
    }
  }  // ::pfasst::encap
} // ::pfasst
