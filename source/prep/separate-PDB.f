C======================================================================C
C                                                                      C
C         This program separate to  single coordinate PDBs             C
C              from multi coordinate PDB data made by VMD soft         C
C                                                                      C
C                          Programed by J. Koseki  @ Nagoya Univ.      C
C                                       (enhanced from split.f)        C
C                                        Version 1.0   15.06.08        C
C                                        Version 1.1   19.08.04        C
C                                                                      C
C======================================================================C

      Implicit None

      Integer         Num,        Keta,
     *                i_check1,   i_check2,   i_check3
      Character*5     C_Num
      Character*70    N_Input,    Title,      Line,      Clear_Data

      Character*20    N_Out

      Write(*, *) "Plese, set your input name. (Asum-PDB-type)"
      Read (*, Fmt = "(A70)") N_Input

      N_Input = Trim(Adjustl(N_Input))

      Open (10, file = N_Input, status = 'old', Iostat = i_check1)
      If (i_check1 /=  0) Then
       Write (*, *) "Input pdb data does not exist."
       Stop
      End if

      Read (10, Fmt = "(A70)") Title
      


      Num = 1

      Do
       Keta = Int(log10(dble(Num))) + 1

       Write(C_Num, Fmt = "(I5)") Num

       Select Case(Keta)
        Case(1)
         C_Num = "0000"//Trim(Adjustl(C_Num))
        Case(2)
         C_Num = "000"//Trim(Adjustl(C_Num))
        Case(3)
         C_Num = "00"//Trim(Adjustl(C_Num))
        Case(4)
         C_Num = "0"//Trim(Adjustl(C_Num))
        Case(5)
         Continue
       End select

       N_Out = "Separate_"//C_Num//".pdb"

       Open (20, file = N_Out, status = 'new', Iostat = i_check2)
       If (i_check2 /=  0) Then
        Write (*, *) "Wirting Error!! in ", N_Out
        Stop
       End if


       Write(* , *)             "During write ", N_Out
       Write(20, Fmt = "(A70)") Title

       Do
        Read (10, Fmt = "(A70)", Iostat = i_check3) Line

        If (i_check3 /= 0) Then
         Close(10)
         Close(20)
         Clear_Data = "rm "//Trim(Adjustl(N_Out))
         Call system(Clear_Data)
         Stop
        End if

        Write(20, Fmt = "(A70)") Line
        If (Index(Line, "END") /= 0) Exit
       End do


       Num = Num + 1


       Close(20)

      End do

      Close(10)

      End
