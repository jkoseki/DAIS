#!/usr/bin/bash

function usage {
        cat <<EOM

Usage: $(basename "$0") -i [PDB Name]
        -i      Input PDB       (ex. Protein.pdb)

EOM

        exit 2
}

PATH=${PATH}:/opt/schrodinger2021-2/:/opt/schrodinger2021-2/utilities
export PATH


while getopts ":i:h" optKey;
do
        case "$optKey" in
                i)
                        pdbname=${OPTARG%.*}
                        intname=$pdbname'-int.pdb'
                        ppaname=$pdbname'-out.pdb'
                        maename=$pdbname'-out.mae'
                        logname=$pdbname'-int.log'
                        grep -v '   H ' $OPTARG | grep -v 'CONECT' > $intname
                        prepwizard -rehtreat -c -disulfides -noepik -s -propka_pH 7.0 -noimpref $intname $ppaname
                        ppwfrag=0
                        while [ "$ppwfrag" != "1" ]
                        do
                                sleep   3
                                ppwfrag=$(ls $ppaname | wc -l)
                        done
                        prepwizard -rehtreat -c -disulfides -noepik -s -propka_pH 7.0 -noimpref $intname $maename
                        ppwfrag=0
                        while [ $ppwfrag != 1 ]
                        do
                                sleep   3
                                ppwfrag=$(ls $maename | wc -l)
                        done

                        grep -v '   H ' $ppaname | grep -v 'CONECT' > $intname
                        mv $intname $ppaname
                        sed -i "s/NMA/NME/g" $ppaname
                        sed -i "s/CA  NME/CH3 NME/g" $ppaname

                        cyxname=$pdbname'-cyx.dat'
                        hiename=$pdbname'-hie.dat'
                        hipname=$pdbname'-hip.dat'

                        grep 'CYX' $maename > $cyxname
                        grep 'HIE' $maename > $hiename
                        grep 'HIP' $maename > $hipname

                        rm $maename $logname

                        sed -i -e 's/ " ".*//' $cyxname
                        sed -i -e 's/ " ".*//' $hiename
                        sed -i -e 's/ " ".*//' $hipname

                        sed -i -e 's/.*\s\([0-9]*\).*/\1/' $cyxname
                        sed -i -e 's/.*\s\([0-9]*\).*/\1/' $hiename
                        sed -i -e 's/.*\s\([0-9]*\).*/\1/' $hipname

                        unumcyx='UN-'$cyxname
                        unumhie='UN-'$hiename
                        unumhip='UN-'$hipname

                        sort $cyxname | uniq -d > $unumcyx
                        sort $hiename | uniq -d > $unumhie
                        sort $hipname | uniq -d > $unumhip

                        rm $cyxname $hiename $hipname


                        sed -i "s/HIS/HID/g" $ppaname

                        while read line
                        do
                                frpt=0
                                case ${#line} in
                                        0)
                                                break
                                                frpt=1
                                                ;;
                                        1)
                                                tnum='   '$line

                                                ;;
                                        2)
                                                tnum='  '$line
                                                ;;
                                        3)
                                                tnum=' '$line
                                                ;;
                                        4)
                                                tnum=$line
                                                ;;
                                        *)
                                                tnum='****'
                                                ;;
                                esac

                                if [ $frpt = 1 ] ; then
                                        break
                                fi

                                beft='CYS A'$tnum
                                aftt='CYX A'$tnum
                                sed -i "s/$beft/$aftt/g" $ppaname

                        done < $unumcyx



                        while read line
                        do
                                frpt=0
                                case ${#line} in
                                        0)
                                                break
                                                frpt=1
                                                ;;
                                        1)
                                                tnum='   '$line

                                                ;;
                                        2)
                                                tnum='  '$line
                                                ;;
                                        3)
                                                tnum=' '$line
                                                ;;
                                        4)
                                                tnum=$line
                                                ;;
                                        *)
                                                tnum='****'
                                                ;;
                                esac

                                if [ $frpt = 1 ] ; then
                                        break
                                fi

                                beft='HID A'$tnum
                                aftt='HIE A'$tnum
                                sed -i "s/$beft/$aftt/g" $ppaname

                        done < $unumhie


                        while read line
                        do
                                frpt=0
                                case ${#line} in
                                        0)
                                                break
                                                frpt=1
                                                ;;
                                        1)
                                                tnum='   '$line

                                                ;;
                                        2)
                                                tnum='  '$line
                                                ;;
                                        3)
                                                tnum=' '$line
                                                ;;
                                        4)
                                                tnum=$line
                                                ;;
                                        *)
                                                tnum='****'
                                                ;;
                                esac

                                if [ $frpt = 1 ] ; then
                                        break
                                fi

                                beft='HID A'$tnum
                                aftt='HIP A'$tnum
                                sed -i "s/$beft/$aftt/g" $ppaname

                        done < $unumhip

                        rm $unumcyx $unumhie $unumhip



                        prmname=$pdbname'.prmtop'
                        crdname=$pdbname'.inpcrd'
                        trpscrp=$pdbname'.scrp'
                        cat <<- EOF > $trpscrp
                                source oldff/leaprc.ff99SB
                                AFP=loadpdb $ppaname
                                addions AFP Na+ 0
                                addions AFP Cl- 0
                                solvatebox AFP TIP3PBOX 12.0
                                saveamberparm AFP $prmname $crdname
                                quit
                        EOF

                        tleap -s -f $trpscrp

                        rm leap.log $ppaname $trpscrp



                        grotop=$pdbname'_GMX.top'
                        grogro=$pdbname'_GMX.gro'

                        acpype.py -p $prmname -x $crdname

                        sleep 5

                        if [ ! -d Amber-input ]; then
                                mkdir Amber-input
                        fi

                        mv $prmname ./Amber-input/.
                        mv $crdname ./Amber-input/.


                        if [ ! -d Gromacs-input ]; then
                                mkdir Gromacs-input
                        fi

                        mv $grotop ./Gromacs-input/.
                        mv $grogro ./Gromacs-input/.

                        rm em.mdp md.mdp

                        ;;
                '-h'|'--help'|* )
                        usage
                        ;;
        esac
done



