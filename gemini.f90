program gemini
   implicit none
   integer :: a(10) = (/5, 1, 5, 8, 5, 10, 6, 2, 5, 7/)
   integer :: a_temporary(10)

   !一時変数に計算結果を記録
   a_temporary = a ** 2 + 3 * a

   !一時的な配列から元の配列に一つ要素をずらしながら入れていく
   a = cshift(a_temporary, shift=-1)

   !出力
   write(*, *) a
   stop
end program gemini