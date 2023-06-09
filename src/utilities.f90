! Reference https://gist.github.com/t-nissie/479f0f16966925fa29ea
recursive subroutine quicksort_int(a, first, last)
   implicit none
   integer  a(*), x, t
   integer first, last
   integer i, j

   x = a( (first+last) / 2 )
   i = first
   j = last
   do
      do while (a(i) < x)
         i=i+1
      end do
      do while (x < a(j))
         j=j-1
      end do
      if (i >= j) exit
      t = a(i);  a(i) = a(j);  a(j) = t
      i=i+1
      j=j-1
   end do
   if (first < i-1) call quicksort_int(a, first, i-1)
   if (j+1 < last)  call quicksort_int(a, j+1, last)
end subroutine quicksort_int