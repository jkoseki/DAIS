#!/usr/bin/bash

function usage {
        cat <<EOM

Usage: $(basename "$0") -i [mdcrd Name header]
        -i      Input mdcrd     (ex. XXX only for XXX.mdcrd)

EOM
        exit 2
}


while getopts ":i:h" optKey;
do
        case "$optKey" in
                i)
                        headname=${OPTARG%}
                        pmtpname=$headname'.prmtop'
                        mdcdname=$headname'.mdcrd'
                        opdbname=$headname'.pdb'
                        notwname='not_water-'$opdbname
                        ipdbname='Int-'$opdbname

                        cpptraj << EOF
                                parm $pmtpname
                                trajin $mdcdname 1 last 1
                                center :1-420
                                trajout $opdbname
EOF

                        grep -v " WAT"   $opdbname > $notwname
                        grep -v "  C   " $notwname > $opdbname
                        grep -v "  O "   $opdbname > $ipdbname
                        grep -v "  N "   $ipdbname > $opdbname
                        grep -v "  H "   $opdbname > $ipdbname
                        grep -v "  HH"   $ipdbname > $opdbname
                        grep -v "  HA"   $opdbname > $ipdbname
                        grep -v "  HB"   $ipdbname > $opdbname
                        grep -v "  HD"   $opdbname > $ipdbname
                        grep -v "  CH"   $ipdbname > $opdbname
                        grep -v "  CB"   $opdbname > $ipdbname
                        grep -v "  CE"   $ipdbname > $opdbname
                        grep -v "  CZ"   $opdbname > $ipdbname
                        grep -v "  CD"   $ipdbname > $opdbname
                        grep -v "  CG"   $opdbname > $ipdbname
                        grep -v "  SD"   $ipdbname > $opdbname
                        grep -v "  SG"   $opdbname > $ipdbname
                        grep -v " Cl-"   $ipdbname > $opdbname
                        grep -v " Na+"   $opdbname > $ipdbname

                        mv $ipdbname $opdbname

                        mkdir $headname
                        mkdir $headname/CA-data
                        mkdir $headname/Location_detect
                        mkdir $headname/Polar-loc_2
                        mkdir $headname/Pu-loc_2

                        sed -i 's/ENDMDL/END/' $notwname
                        sed -i 's/ENDMDL/END/' $opdbname

                        sed -i 's/END   //' $notwname
                        sed -i 's/END   //' $opdbname

                        mv $notwname $headname/.
                        mv $opdbname $headname/CA-data/.


                        cd $headname/;separate-pdb2.exe << EOF
$notwname
EOF

                        rm $notwname

                        cd CA-data/;separate-pdb2.exe << EOF
$opdbname
EOF

                        rm $opdbname

                        cd ../../

                        ;;
                '-h'|'--help'|* )
                        usage
                        ;;
        esac
done

