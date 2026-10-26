module sha256_mod
   use iso_fortran_env, only: int8, int32, int64
   implicit none
   private
   public :: sha256_file
   integer(int64), parameter :: mask32 = int(z'FFFFFFFF', int64)
   integer(int64), parameter :: k(64) = [ &
      int(z'428A2F98',int64),int(z'71374491',int64),int(z'B5C0FBCF',int64),int(z'E9B5DBA5',int64), &
      int(z'3956C25B',int64),int(z'59F111F1',int64),int(z'923F82A4',int64),int(z'AB1C5ED5',int64), &
      int(z'D807AA98',int64),int(z'12835B01',int64),int(z'243185BE',int64),int(z'550C7DC3',int64), &
      int(z'72BE5D74',int64),int(z'80DEB1FE',int64),int(z'9BDC06A7',int64),int(z'C19BF174',int64), &
      int(z'E49B69C1',int64),int(z'EFBE4786',int64),int(z'0FC19DC6',int64),int(z'240CA1CC',int64), &
      int(z'2DE92C6F',int64),int(z'4A7484AA',int64),int(z'5CB0A9DC',int64),int(z'76F988DA',int64), &
      int(z'983E5152',int64),int(z'A831C66D',int64),int(z'B00327C8',int64),int(z'BF597FC7',int64), &
      int(z'C6E00BF3',int64),int(z'D5A79147',int64),int(z'06CA6351',int64),int(z'14292967',int64), &
      int(z'27B70A85',int64),int(z'2E1B2138',int64),int(z'4D2C6DFC',int64),int(z'53380D13',int64), &
      int(z'650A7354',int64),int(z'766A0ABB',int64),int(z'81C2C92E',int64),int(z'92722C85',int64), &
      int(z'A2BFE8A1',int64),int(z'A81A664B',int64),int(z'C24B8B70',int64),int(z'C76C51A3',int64), &
      int(z'D192E819',int64),int(z'D6990624',int64),int(z'F40E3585',int64),int(z'106AA070',int64), &
      int(z'19A4C116',int64),int(z'1E376C08',int64),int(z'2748774C',int64),int(z'34B0BCB5',int64), &
      int(z'391C0CB3',int64),int(z'4ED8AA4A',int64),int(z'5B9CCA4F',int64),int(z'682E6FF3',int64), &
      int(z'748F82EE',int64),int(z'78A5636F',int64),int(z'84C87814',int64),int(z'8CC70208',int64), &
      int(z'90BEFFFA',int64),int(z'A4506CEB',int64),int(z'BEF9A3F7',int64),int(z'C67178F2',int64) ]
contains
   pure function rotr(x, n) result(y)
      integer(int64), intent(in) :: x
      integer, intent(in) :: n
      integer(int64) :: y
      y = iand(ior(ishft(x, -n), ishft(x, 32-n)), mask32)
   end function rotr

   function sha256_file(path) result(hex)
      character(len=*), intent(in) :: path
      character(len=64) :: hex
      integer(int8), allocatable :: bytes(:)
      integer(int64) :: h(8), w(64), a,b,c,d,e,f,g,hh,t1,t2,s0,s1,ch,maj
      integer(int64) :: n, padded, block, pos, bits, x1,x2,x3,x4
      integer :: unit, ios, i, j, idx

      hex = ''
      inquire(file=path, size=n, iostat=ios)
      if (ios /= 0 .or. n < 0) return
      padded = ((n + 9 + 63) / 64) * 64
      allocate(bytes(padded), stat=ios)
      if (ios /= 0) return
      bytes = 0_int8
      open(newunit=unit, file=path, access='stream', form='unformatted', status='old', action='read', iostat=ios)
      if (ios /= 0) return
      if (n > 0) read(unit, iostat=ios) bytes(1:n)
      close(unit)
      if (ios /= 0) return
      bytes(n+1) = int(z'80', int8)
      bits = n*8
      do i = 0, 7
         bytes(padded-i) = int(iand(ishft(bits, -8*i), 255_int64), int8)
      end do

      h = [int(z'6A09E667',int64),int(z'BB67AE85',int64),int(z'3C6EF372',int64),int(z'A54FF53A',int64), &
           int(z'510E527F',int64),int(z'9B05688C',int64),int(z'1F83D9AB',int64),int(z'5BE0CD19',int64)]
      do block = 0, padded/64-1
         pos = block*64
         do i = 1, 16
            idx = int(pos + (i-1)*4 + 1)
            x1 = iand(int(bytes(idx),int64),255_int64)
            x2 = iand(int(bytes(idx+1),int64),255_int64)
            x3 = iand(int(bytes(idx+2),int64),255_int64)
            x4 = iand(int(bytes(idx+3),int64),255_int64)
            w(i) = iand(ishft(x1,24)+ishft(x2,16)+ishft(x3,8)+x4,mask32)
         end do
         do i = 17, 64
            s0 = ieor(ieor(rotr(w(i-15),7),rotr(w(i-15),18)),ishft(w(i-15),-3))
            s1 = ieor(ieor(rotr(w(i-2),17),rotr(w(i-2),19)),ishft(w(i-2),-10))
            w(i) = iand(w(i-16)+s0+w(i-7)+s1,mask32)
         end do
         a=h(1); b=h(2); c=h(3); d=h(4); e=h(5); f=h(6); g=h(7); hh=h(8)
         do i = 1, 64
            s1 = ieor(ieor(rotr(e,6),rotr(e,11)),rotr(e,25))
            ch = ieor(iand(e,f),iand(ieor(e,mask32),g))
            t1 = iand(hh+s1+ch+k(i)+w(i),mask32)
            s0 = ieor(ieor(rotr(a,2),rotr(a,13)),rotr(a,22))
            maj = ieor(ieor(iand(a,b),iand(a,c)),iand(b,c))
            t2 = iand(s0+maj,mask32)
            hh=g; g=f; f=e; e=iand(d+t1,mask32); d=c; c=b; b=a; a=iand(t1+t2,mask32)
         end do
         h = [iand(h(1)+a,mask32),iand(h(2)+b,mask32),iand(h(3)+c,mask32),iand(h(4)+d,mask32), &
               iand(h(5)+e,mask32),iand(h(6)+f,mask32),iand(h(7)+g,mask32),iand(h(8)+hh,mask32)]
      end do
      do i = 1, 8
         write(hex((i-1)*8+1:i*8),'(z8.8)') int(h(i), int32)
      end do
      do i = 1, 64
         if (hex(i:i) >= 'A' .and. hex(i:i) <= 'F') hex(i:i) = achar(iachar(hex(i:i))+32)
      end do
      deallocate(bytes)
   end function sha256_file
end module sha256_mod
