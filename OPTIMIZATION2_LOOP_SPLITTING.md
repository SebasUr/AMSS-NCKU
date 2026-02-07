# AMSS-NCKU Loop Splitting Optimization

## Summary

**Optimization:** Split interior/boundary loops in `diff_new.f90` to remove conditionals from hot loops.

**Results:**
- Before (commit c290a95): 37,185 seconds
- After (Job 1178): **36,924 seconds**
- **Improvement: 261 seconds faster (0.7%)**

**Cumulative speedup from baseline:** ~1.05x + 0.7% = **~1.06x**

---

## What Was Changed

### File: `AMSS_NCKU_source/diff_new.f90`

**Routines optimized:** `fderivs` and `fdderivs` (4th order section, `ghost_width == 3`)

**Before:**
```fortran
do k=1,ex(3)-1
do j=1,ex(2)-1
do i=1,ex(1)-1
   if(i+2 <= imax .and. i-2 >= imin .and. &
      j+2 <= jmax .and. j-2 >= jmin .and. &
      k+2 <= kmax .and. k-2 >= kmin) then
      ! 4th order stencil
   elseif(...) then
      ! 2nd order fallback
   endif
enddo
enddo
enddo
```

**After:**
```fortran
! Interior loop - no conditionals, fully vectorizable
do k=3,ex(3)-3
do j=3,ex(2)-3
do i=3,ex(1)-3
    fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))
    fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))
    fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))
enddo
enddo
enddo

! Boundary loop - only ~5% of points, with conditionals
do k=1,ex(3)-1
do j=1,ex(2)-1
do i=1,ex(1)-1
  if(i < 3 .or. i > ex(1)-3 .or. j < 3 .or. j > ex(2)-3 .or. k < 3 .or. k > ex(3)-3) then
    ! original conditional logic for boundaries
  endif
enddo
enddo
enddo
```

---

## Why This Works

1. **Removes branch misprediction:** The interior loop has no `if` statements, allowing CPU to execute without stalls.

2. **Enables vectorization:** Without conditionals, the compiler can generate SIMD instructions (AVX2/AVX-512) for the inner loop.

3. **Better cache utilization:** Predictable memory access pattern in interior allows hardware prefetcher to work efficiently.

4. **~95% of points are interior:** Only boundary points (2-3 layers) need the conditional logic.

---

## Compliance with Competition Rules

✅ **VALID**: This is a code implementation optimization only.

- ❌ **NO changes** to numerical method (still 4th order finite differences)
- ❌ **NO changes** to input parameters
- ❌ **NO changes** to results (mathematically equivalent)
- ✅ **ONLY changes**: Loop structure for better vectorization

---

## Files Modified

- `AMSS_NCKU_source/diff_new.f90` (loop splitting in `fderivs` and `fdderivs`)

---

**Date:** February 7, 2026
**Job ID:** 1178
**Result:** 0.7% additional speedup (37,185s → 36,924s)
