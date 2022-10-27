!**********************************************************************!
!                                                                      !
!   This program is many mutation type protein script maker            !
!                                  with protein preparation wizard.    !
!                                                                      !
!                                     Programed by J. Koseki           !
!                                     V1.0  21.10.28 @ Nagoya Univ.    !
!                                                                      !
!**********************************************************************!

      Implicit None

      Integer    i,            j,            k,           l,           &
                 m,            n,            o,                        &
                 Raw_Num,      Int_Step,     All_Pattern, Sep,         &
                 Pre_Sep,      Target_reg,                             &
                 icheck_1,     icheck_2,     icheck_3

      Integer,                 Allocatable,  Dimension(:) :: Index_Num,&
                                                             MPLen_Num


      Character*1,             Allocatable,  Dimension(:, :) :: MP_Mat

      Character*30,            Allocatable,  Dimension(:) :: Mut_Patt

      Character*3              Det_TEXT,     Reg_Name^M
      Character*100            Input_Name,   Output_Name,  Target_Name
      Character*200            PATH



      !  ===  Load Mutation infomations
      Write(*, Fmt = "('============================================')")
      Write(*, Fmt = "(' Start program to make Script ')")
      Write(*, Fmt = "('   for make Mutation body for Schrodinger')")
      Write(*, Fmt = "('    Protein preparation Wizard in Maestro.')")
      Write(*, Fmt = "('============================================')")
      Write(*, Fmt = "('')")
      Write(*, Fmt = "('--- Load Mutation information  ---')")

      Input_Name = "Mut-info"
      Open(10, file = Input_Name, status = 'old', Iostat = icheck_1)

      If (icheck_1 /= 0 ) Then
          Write(*, Fmt = "('Plese, check your Mut-info file.')")
          Write(*, Fmt = "('Mut-info file could not be found.')")
          Stop
      End if


      Index_Num(:) = 0
      Int_Step     = 0

      Read(10, *)

      Do
          Read(10, Fmt ="(I)", Iostat = icheck_2) Raw_Num
          If (icheck_2 /= 0) Exit^M
          If (Raw_Num == 0) Exit ^M
          Int_Step = Int_Step + 1
      End do

      Close(10)
      Allocate(Index_Num(Int_Step), Mut_Patt(Int_Step),                &
               MPLen_Num(Int_Step))
      Open(10, file = Input_Name)


      Read(10, *)
      Do i = 1, Int_Step
          Read(10, Fmt = "(I)") Index_Num(i)
      End do
      Read(10, *)
      Read(10, *)
      Do j = 1, Int_Step
          Read(10, Fmt = "(A)") Mut_Patt(j)
          MPLen_Num(j) = Len_Trim(Mut_Patt(j))
      End do

      Close(10)




      ! Make Script pattern of mutation

      All_Pattern = 1

      Do k = 1, Int_Step
          All_Pattern = All_Pattern * MPLen_Num(k)
      End do

      Allocate(MP_Mat(All_Pattern, Int_Step))


      Pre_Sep = All_Pattern

      Do l = 1, Int_Step
          Sep = Pre_Sep / MPLen_Num(l)

          Do m = 1, All_Pattern
              Target_reg = 1 + ((m - 1) / Sep) - MPLen_Num(l) * ((m - 1) / Sep / MPLen_Num(l))
              MP_Mat(m, l) = Mut_Patt(l)(Target_reg:Target_reg)
          End do

          Pre_Sep = Sep

      End do


      ! Write script
      Input_Name = "Mut-Script.cmd"
      Open(20, file = Input_Name, status = 'new', Iostat = icheck_1)

      If (icheck_1 /= 0 ) Then
          Write(*, Fmt = "('Plese, check your Mut-info file.')")
          Write(*, Fmt = "('Mut-info file could not be found.')")
          Stop
      End if


      Do n = 1, All_Pattern
          PATH = '"./'

          Do o = 1, Int_Step
              Write(20, Fmt = "('workspaceselectionreplace fillres atom.num ', I4)") Index_Num(o)
              Write(20, Fmt = "('prefer  wsselectpickstate=residue')")

              Select Case(MP_Mat(n, o))
                  Case ('A')
                      Reg_Name = 'ALA'
                  Case ('R')
                      Reg_Name = 'ARG'
                  Case ('N')
                      Reg_Name = 'ASN'
                  Case ('D')
                      Reg_Name = 'ASP'
                  Case ('C')
                      Reg_Name = 'CYS'
                  Case ('Q')
                      Reg_Name = 'GLN'
                  Case ('E')
                      Reg_Name = 'GLU'
                  Case ('G')
                      Reg_Name = 'GLY'
                  Case ('H')
                      Reg_Name = 'HIS'
                  Case ('I')
                      Reg_Name = 'ILE'
                  Case ('L')
                      Reg_Name = 'LEU'
                  Case ('K')
                      Reg_Name = 'LYS'
                  Case ('M')
                      Reg_Name = 'MET'
                  Case ('F')
                      Reg_Name = 'PHE'
                  Case ('P')
                      Reg_Name = 'PRO'
                  Case ('S')
                      Reg_Name = 'SER'
                  Case ('T')
                      Reg_Name = 'THR'
                  Case ('W')
                      Reg_Name = 'TRP'
                  Case ('Y')
                      Reg_Name = 'TYR'
                  Case ('V')
                      Reg_Name = 'VAL'
              End select

              Write(20, Fmt = "('fragment peptide ', A3)") Reg_Name
              Write(20, Fmt = "('mutate workspace_selection and protein')")

              PATH = Trim(PATH)//MP_Mat(n, o)
          End do

          PATH = Trim(PATH)//'.pdb"'

          Write(20, Fmt = "('')")
          Write(20, Fmt = "('highlighttextunselectall')")
          Write(20, Fmt = "('minimize atom.selected ')")
          Write(20, Fmt = "('entryexport  format=pdb')")
          Write(20, Fmt = "('entryexport  format=xyz')")
          Write(20, Fmt = "('entryexport  format=pdb')")
          Write(20, Fmt = "('entryexport ', A15)") Adjustl(PATH)
          Write(20, Fmt = "('')")

      End do

      Close(20)

      End
