#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --title "AUTOMATIC WAVE FUNCTION GENERATER IN VASP_CHGCAR FORMAT" --colors --infobox "\n\Zb\Z1
	=======================================================================================\n
	 THIS PROGRAM WILL GENERATE SAPERATE EIGEN VECTOR FILES FROM THE COMMON OUTPUT FILES \n 
                                     IN A GIVEN RANGE \n
	=======================================================================================\n \Z0 
	\n\n          GIVE THE FOLLOWING PARAMETERS SEPERATED BY A <SPACE> \n\n
                         WHICH SPIN EIGEN VECTORS YOU WANT (UP/DN)?\n
                         GIVE THE LOWEST VECTOR YOUR WANT::\n
                         GIVE THE HIGHEST VECTOR YOUR WANT::\n\n\n
                       \Z1\Zu TO GO BACK TO MAIN MENU TYPE q" 20 100
read spin vecmin vecmax
case "$spin"
in
q) inf=2;dialog --backtitle "PROGRAM BY SANDIP DE   DATE : 19.02.09 " --colors --timeout 1 --msgbox " \n \Zb\Z1\Zu GOING BACK TO MAIN MENU " 10 40;;
*)

for vec in `seq ${vecmin} 1 ${vecmax}`;
do
eigenvecfile='1st6_eigen_vec'
if [ ${vec} -gt 6 ]
then
  eigenvecfile='7to12_eigen_vec'
  if [ ${vec} -gt 12 ]
	then
	eigenvecfile='13to18_eigen_vec'
	if [ ${vec} -gt 18 ]
	then
	  eigenvecfile='19to24_eigen_vec'
	  if [ ${vec} -gt 24 ]
	  then
	    eigenvecfile='25to30_eigen_vec' 
	   fi
	fi	
  fi
fi 
cat >awk.awk<< EOF
	BEGIN{
	  col=vec-6*int(vec/6)
	  if ( col == 0 ) col=6
	}
	{i=col}{print $ i}
EOF
cat >eigen_vec_spin${spin}${vec}.VASP_CHGCAR <vmdformat 
awk -v vec=${vec} -f awk.awk <${eigenvecfile}${spin}.dat >>eigen_vec_spin${spin}${vec}.VASP_CHGCAR
if [ ${vec} -eq 1 ]
then
	printf "Your  ${vec} st spin${spin} eigen-vector is stored in file:eigen_vec_spin${spin}${vec}.VASP_CHGCAR from ${eigenvecfile}${spin}.dat\n"
fi
if [ ${vec} -eq 2 ]
then
	printf "Your  ${vec} nd spin${spin}eigen-vector is stored in file:eigen_vec_spin${spin}${vec}.VASP_CHGCAR from ${eigenvecfile}${spin}.dat\n"
fi
if [ ${vec} -eq 3 ]
then
	printf "Your  ${vec} rd spin${spin} eigen-vector is stored in file:eigen_vec_spin${spin}${vec}.VASP_CHGCAR from ${eigenvecfile}${spin}.dat\n"
fi
if [ ${vec} -gt 3 ]
then
	printf "Your  ${vec} th spin${spin}eigen-vector is stored in file:eigen_vec_spin${spin}${vec}.VASP_CHGCAR from ${eigenvecfile}${spin}.dat\n"

fi
done
rm -f awk.awk
dialog --msgbox " OUTPUT GENERATED SUCCESSFULLY " 8 40 ;;
esac
done
