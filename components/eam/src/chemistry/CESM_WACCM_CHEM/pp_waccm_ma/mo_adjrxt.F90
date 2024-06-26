      module mo_adjrxt
      private
      public :: adjrxt
      contains
      subroutine adjrxt( rate, inv, m, ncol, nlev )
      use shr_kind_mod, only : r8 => shr_kind_r8
      use chem_mods, only : nfs, rxntot
      implicit none
!--------------------------------------------------------------------
! ... dummy arguments
!--------------------------------------------------------------------
      integer, intent(in) :: ncol, nlev
      real(r8), intent(in) :: inv(ncol,nlev,nfs)
      real(r8), intent(in) :: m(ncol,nlev)
      real(r8), intent(inout) :: rate(ncol,nlev,rxntot)
      rate(:,:, 91) = rate(:,:, 91) * inv(:,:, 2)
      rate(:,:, 95) = rate(:,:, 95) * inv(:,:, 2)
      rate(:,:, 99) = rate(:,:, 99) * inv(:,:, 2)
      rate(:,:, 104) = rate(:,:, 104) * inv(:,:, 1)
      rate(:,:, 105) = rate(:,:, 105) * inv(:,:, 1)
      rate(:,:, 111) = rate(:,:, 111) * inv(:,:, 1)
      rate(:,:, 121) = rate(:,:, 121) * inv(:,:, 1)
      rate(:,:, 133) = rate(:,:, 133) * inv(:,:, 1)
      rate(:,:, 141) = rate(:,:, 141) * inv(:,:, 1)
      rate(:,:, 144) = rate(:,:, 144) * inv(:,:, 1)
      rate(:,:, 145) = rate(:,:, 145) * inv(:,:, 1)
      rate(:,:, 146) = rate(:,:, 146) * inv(:,:, 1)
      rate(:,:, 148) = rate(:,:, 148) * inv(:,:, 1)
      rate(:,:, 149) = rate(:,:, 149) * inv(:,:, 1)
      rate(:,:, 164) = rate(:,:, 164) * inv(:,:, 1)
      rate(:,:, 184) = rate(:,:, 184) * inv(:,:, 1)
      rate(:,:, 185) = rate(:,:, 185) * inv(:,:, 1)
      rate(:,:, 195) = rate(:,:, 195) * inv(:,:, 1)
      rate(:,:, 237) = rate(:,:, 237) * inv(:,:, 1)
      rate(:,:, 276) = rate(:,:, 276) * inv(:,:, 2)
      rate(:,:, 279) = rate(:,:, 279) * inv(:,:, 2)
      rate(:,:, 283) = rate(:,:, 283) * inv(:,:, 2)
      rate(:,:, 288) = rate(:,:, 288) * inv(:,:, 2)
      rate(:,:, 289) = rate(:,:, 289) * inv(:,:, 2)
      rate(:,:, 89) = rate(:,:, 89) * m(:,:)
      rate(:,:, 90) = rate(:,:, 90) * m(:,:)
      rate(:,:, 92) = rate(:,:, 92) * m(:,:)
      rate(:,:, 93) = rate(:,:, 93) * m(:,:)
      rate(:,:, 94) = rate(:,:, 94) * m(:,:)
      rate(:,:, 96) = rate(:,:, 96) * m(:,:)
      rate(:,:, 97) = rate(:,:, 97) * m(:,:)
      rate(:,:, 98) = rate(:,:, 98) * m(:,:)
      rate(:,:, 100) = rate(:,:, 100) * m(:,:)
      rate(:,:, 101) = rate(:,:, 101) * m(:,:)
      rate(:,:, 102) = rate(:,:, 102) * m(:,:)
      rate(:,:, 103) = rate(:,:, 103) * m(:,:)
      rate(:,:, 104) = rate(:,:, 104) * m(:,:)
      rate(:,:, 105) = rate(:,:, 105) * m(:,:)
      rate(:,:, 106) = rate(:,:, 106) * m(:,:)
      rate(:,:, 107) = rate(:,:, 107) * m(:,:)
      rate(:,:, 108) = rate(:,:, 108) * m(:,:)
      rate(:,:, 109) = rate(:,:, 109) * m(:,:)
      rate(:,:, 110) = rate(:,:, 110) * m(:,:)
      rate(:,:, 111) = rate(:,:, 111) * m(:,:)
      rate(:,:, 112) = rate(:,:, 112) * m(:,:)
      rate(:,:, 113) = rate(:,:, 113) * m(:,:)
      rate(:,:, 114) = rate(:,:, 114) * m(:,:)
      rate(:,:, 115) = rate(:,:, 115) * m(:,:)
      rate(:,:, 116) = rate(:,:, 116) * m(:,:)
      rate(:,:, 117) = rate(:,:, 117) * m(:,:)
      rate(:,:, 118) = rate(:,:, 118) * m(:,:)
      rate(:,:, 119) = rate(:,:, 119) * m(:,:)
      rate(:,:, 120) = rate(:,:, 120) * m(:,:)
      rate(:,:, 121) = rate(:,:, 121) * m(:,:)
      rate(:,:, 122) = rate(:,:, 122) * m(:,:)
      rate(:,:, 123) = rate(:,:, 123) * m(:,:)
      rate(:,:, 124) = rate(:,:, 124) * m(:,:)
      rate(:,:, 125) = rate(:,:, 125) * m(:,:)
      rate(:,:, 126) = rate(:,:, 126) * m(:,:)
      rate(:,:, 127) = rate(:,:, 127) * m(:,:)
      rate(:,:, 128) = rate(:,:, 128) * m(:,:)
      rate(:,:, 129) = rate(:,:, 129) * m(:,:)
      rate(:,:, 130) = rate(:,:, 130) * m(:,:)
      rate(:,:, 131) = rate(:,:, 131) * m(:,:)
      rate(:,:, 132) = rate(:,:, 132) * m(:,:)
      rate(:,:, 133) = rate(:,:, 133) * m(:,:)
      rate(:,:, 134) = rate(:,:, 134) * m(:,:)
      rate(:,:, 135) = rate(:,:, 135) * m(:,:)
      rate(:,:, 136) = rate(:,:, 136) * m(:,:)
      rate(:,:, 137) = rate(:,:, 137) * m(:,:)
      rate(:,:, 138) = rate(:,:, 138) * m(:,:)
      rate(:,:, 139) = rate(:,:, 139) * m(:,:)
      rate(:,:, 140) = rate(:,:, 140) * m(:,:)
      rate(:,:, 141) = rate(:,:, 141) * m(:,:)
      rate(:,:, 142) = rate(:,:, 142) * m(:,:)
      rate(:,:, 143) = rate(:,:, 143) * m(:,:)
      rate(:,:, 144) = rate(:,:, 144) * m(:,:)
      rate(:,:, 145) = rate(:,:, 145) * m(:,:)
      rate(:,:, 146) = rate(:,:, 146) * m(:,:)
      rate(:,:, 147) = rate(:,:, 147) * m(:,:)
      rate(:,:, 150) = rate(:,:, 150) * m(:,:)
      rate(:,:, 151) = rate(:,:, 151) * m(:,:)
      rate(:,:, 152) = rate(:,:, 152) * m(:,:)
      rate(:,:, 153) = rate(:,:, 153) * m(:,:)
      rate(:,:, 154) = rate(:,:, 154) * m(:,:)
      rate(:,:, 155) = rate(:,:, 155) * m(:,:)
      rate(:,:, 156) = rate(:,:, 156) * m(:,:)
      rate(:,:, 157) = rate(:,:, 157) * m(:,:)
      rate(:,:, 158) = rate(:,:, 158) * m(:,:)
      rate(:,:, 159) = rate(:,:, 159) * m(:,:)
      rate(:,:, 160) = rate(:,:, 160) * m(:,:)
      rate(:,:, 161) = rate(:,:, 161) * m(:,:)
      rate(:,:, 162) = rate(:,:, 162) * m(:,:)
      rate(:,:, 163) = rate(:,:, 163) * m(:,:)
      rate(:,:, 164) = rate(:,:, 164) * m(:,:)
      rate(:,:, 165) = rate(:,:, 165) * m(:,:)
      rate(:,:, 166) = rate(:,:, 166) * m(:,:)
      rate(:,:, 167) = rate(:,:, 167) * m(:,:)
      rate(:,:, 168) = rate(:,:, 168) * m(:,:)
      rate(:,:, 169) = rate(:,:, 169) * m(:,:)
      rate(:,:, 170) = rate(:,:, 170) * m(:,:)
      rate(:,:, 171) = rate(:,:, 171) * m(:,:)
      rate(:,:, 172) = rate(:,:, 172) * m(:,:)
      rate(:,:, 173) = rate(:,:, 173) * m(:,:)
      rate(:,:, 174) = rate(:,:, 174) * m(:,:)
      rate(:,:, 175) = rate(:,:, 175) * m(:,:)
      rate(:,:, 176) = rate(:,:, 176) * m(:,:)
      rate(:,:, 177) = rate(:,:, 177) * m(:,:)
      rate(:,:, 178) = rate(:,:, 178) * m(:,:)
      rate(:,:, 179) = rate(:,:, 179) * m(:,:)
      rate(:,:, 180) = rate(:,:, 180) * m(:,:)
      rate(:,:, 181) = rate(:,:, 181) * m(:,:)
      rate(:,:, 182) = rate(:,:, 182) * m(:,:)
      rate(:,:, 183) = rate(:,:, 183) * m(:,:)
      rate(:,:, 184) = rate(:,:, 184) * m(:,:)
      rate(:,:, 186) = rate(:,:, 186) * m(:,:)
      rate(:,:, 187) = rate(:,:, 187) * m(:,:)
      rate(:,:, 188) = rate(:,:, 188) * m(:,:)
      rate(:,:, 189) = rate(:,:, 189) * m(:,:)
      rate(:,:, 190) = rate(:,:, 190) * m(:,:)
      rate(:,:, 191) = rate(:,:, 191) * m(:,:)
      rate(:,:, 192) = rate(:,:, 192) * m(:,:)
      rate(:,:, 193) = rate(:,:, 193) * m(:,:)
      rate(:,:, 194) = rate(:,:, 194) * m(:,:)
      rate(:,:, 195) = rate(:,:, 195) * m(:,:)
      rate(:,:, 196) = rate(:,:, 196) * m(:,:)
      rate(:,:, 197) = rate(:,:, 197) * m(:,:)
      rate(:,:, 198) = rate(:,:, 198) * m(:,:)
      rate(:,:, 199) = rate(:,:, 199) * m(:,:)
      rate(:,:, 200) = rate(:,:, 200) * m(:,:)
      rate(:,:, 201) = rate(:,:, 201) * m(:,:)
      rate(:,:, 202) = rate(:,:, 202) * m(:,:)
      rate(:,:, 203) = rate(:,:, 203) * m(:,:)
      rate(:,:, 204) = rate(:,:, 204) * m(:,:)
      rate(:,:, 205) = rate(:,:, 205) * m(:,:)
      rate(:,:, 206) = rate(:,:, 206) * m(:,:)
      rate(:,:, 207) = rate(:,:, 207) * m(:,:)
      rate(:,:, 208) = rate(:,:, 208) * m(:,:)
      rate(:,:, 209) = rate(:,:, 209) * m(:,:)
      rate(:,:, 210) = rate(:,:, 210) * m(:,:)
      rate(:,:, 211) = rate(:,:, 211) * m(:,:)
      rate(:,:, 212) = rate(:,:, 212) * m(:,:)
      rate(:,:, 213) = rate(:,:, 213) * m(:,:)
      rate(:,:, 214) = rate(:,:, 214) * m(:,:)
      rate(:,:, 215) = rate(:,:, 215) * m(:,:)
      rate(:,:, 216) = rate(:,:, 216) * m(:,:)
      rate(:,:, 217) = rate(:,:, 217) * m(:,:)
      rate(:,:, 218) = rate(:,:, 218) * m(:,:)
      rate(:,:, 219) = rate(:,:, 219) * m(:,:)
      rate(:,:, 220) = rate(:,:, 220) * m(:,:)
      rate(:,:, 221) = rate(:,:, 221) * m(:,:)
      rate(:,:, 222) = rate(:,:, 222) * m(:,:)
      rate(:,:, 223) = rate(:,:, 223) * m(:,:)
      rate(:,:, 224) = rate(:,:, 224) * m(:,:)
      rate(:,:, 225) = rate(:,:, 225) * m(:,:)
      rate(:,:, 226) = rate(:,:, 226) * m(:,:)
      rate(:,:, 227) = rate(:,:, 227) * m(:,:)
      rate(:,:, 228) = rate(:,:, 228) * m(:,:)
      rate(:,:, 229) = rate(:,:, 229) * m(:,:)
      rate(:,:, 230) = rate(:,:, 230) * m(:,:)
      rate(:,:, 231) = rate(:,:, 231) * m(:,:)
      rate(:,:, 232) = rate(:,:, 232) * m(:,:)
      rate(:,:, 233) = rate(:,:, 233) * m(:,:)
      rate(:,:, 234) = rate(:,:, 234) * m(:,:)
      rate(:,:, 235) = rate(:,:, 235) * m(:,:)
      rate(:,:, 236) = rate(:,:, 236) * m(:,:)
      rate(:,:, 237) = rate(:,:, 237) * m(:,:)
      rate(:,:, 238) = rate(:,:, 238) * m(:,:)
      rate(:,:, 239) = rate(:,:, 239) * m(:,:)
      rate(:,:, 240) = rate(:,:, 240) * m(:,:)
      rate(:,:, 241) = rate(:,:, 241) * m(:,:)
      rate(:,:, 247) = rate(:,:, 247) * m(:,:)
      rate(:,:, 252) = rate(:,:, 252) * m(:,:)
      rate(:,:, 253) = rate(:,:, 253) * m(:,:)
      rate(:,:, 254) = rate(:,:, 254) * m(:,:)
      rate(:,:, 257) = rate(:,:, 257) * m(:,:)
      rate(:,:, 258) = rate(:,:, 258) * m(:,:)
      rate(:,:, 259) = rate(:,:, 259) * m(:,:)
      rate(:,:, 262) = rate(:,:, 262) * m(:,:)
      rate(:,:, 266) = rate(:,:, 266) * m(:,:)
      rate(:,:, 267) = rate(:,:, 267) * m(:,:)
      rate(:,:, 268) = rate(:,:, 268) * m(:,:)
      rate(:,:, 269) = rate(:,:, 269) * m(:,:)
      rate(:,:, 270) = rate(:,:, 270) * m(:,:)
      rate(:,:, 271) = rate(:,:, 271) * m(:,:)
      rate(:,:, 272) = rate(:,:, 272) * m(:,:)
      rate(:,:, 273) = rate(:,:, 273) * m(:,:)
      rate(:,:, 274) = rate(:,:, 274) * m(:,:)
      rate(:,:, 275) = rate(:,:, 275) * m(:,:)
      rate(:,:, 277) = rate(:,:, 277) * m(:,:)
      rate(:,:, 278) = rate(:,:, 278) * m(:,:)
      rate(:,:, 280) = rate(:,:, 280) * m(:,:)
      rate(:,:, 281) = rate(:,:, 281) * m(:,:)
      rate(:,:, 282) = rate(:,:, 282) * m(:,:)
      rate(:,:, 284) = rate(:,:, 284) * m(:,:)
      rate(:,:, 285) = rate(:,:, 285) * m(:,:)
      rate(:,:, 286) = rate(:,:, 286) * m(:,:)
      rate(:,:, 287) = rate(:,:, 287) * m(:,:)
      rate(:,:, 290) = rate(:,:, 290) * m(:,:)
      end subroutine adjrxt
      end module mo_adjrxt
