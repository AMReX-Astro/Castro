
      subroutine get_grav_const(Gconst_out)

         use fundamental_constants_module, only: Gconst

         double precision :: Gconst_out

         Gconst_out = Gconst

      end subroutine get_grav_const

