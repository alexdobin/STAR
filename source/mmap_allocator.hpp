#ifndef _MALLOC_ALLOCATOR_H
#define _MALLOC_ALLOCATOR_H 1

// g++ -std=gnu++20

#include <cstddef>
#include <cstdlib>
#include <new>
#include <type_traits>

#include "mmap_malloc.h"

namespace galaxy {

template <typename _Tp>
class mmap_allocator {
   public:
	using value_type = _Tp;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;
	// _GLIBCXX_RESOLVE_LIB_DEFECTS
	// 2103. propagate_on_container_move_assignment
	using propagate_on_container_move_assignment = std::true_type;

	constexpr mmap_allocator() noexcept = default;
	constexpr mmap_allocator(const mmap_allocator&) noexcept = default;

	template <typename _Tp1>
	constexpr mmap_allocator(const mmap_allocator<_Tp1>&) noexcept {};

	~mmap_allocator() noexcept = default;

	// NB: __n is permitted to be 0.  The C++ standard says nothing
	// about what the return value is when __n == 0.
	[[nodiscard]] _Tp* allocate(size_type __n, const void* = nullptr) {
		// _GLIBCXX_RESOLVE_LIB_DEFECTS
		// 3308. std::allocator<void>().allocate(n)
		static_assert(sizeof(_Tp) != 0, "cannot allocate incomplete types");

		if (__n > this->_M_max_size()) [[unlikely]] {
			// _GLIBCXX_RESOLVE_LIB_DEFECTS
			// 3190. allocator::allocate sometimes returns too little storage
			if (__n > (std::size_t(-1) / sizeof(_Tp))) throw std::bad_array_new_length();
			throw std::bad_alloc();
		}
		_Tp* __ret = static_cast<_Tp*>(mmap_malloc(__n * sizeof(_Tp)));
		if (!__ret) throw std::bad_alloc();
		return __ret;
	}

	// __p is not permitted to be a null pointer.
	void deallocate(_Tp* __p, size_type) { mmap_free(static_cast<void*>(__p)); }

	template <typename _Up>
	friend constexpr bool operator==(const mmap_allocator&, const mmap_allocator<_Up>&) noexcept {
		return true;
	}

   private:
	constexpr size_type _M_max_size() const noexcept {
#if __PTRDIFF_MAX__ < __SIZE_MAX__
		return std::size_t(__PTRDIFF_MAX__) / sizeof(_Tp);
#else
		return std::size_t(-1) / sizeof(_Tp);
#endif
	}
};

}  // namespace galaxy

#endif
