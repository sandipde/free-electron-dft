#!/bin/sh
inf=1
while [ "$inf" = 1 ]
do
zenity  --question --title " 3D_DFT WAVE FUNCTION SEPERATER " --text  "\n 
   =======================================================================================\n
	 THIS PROGRAM WILL GENERATE SAPERATE EIGEN VECTOR FILES FROM THE COMMON OUTPUT FILES \n 
                                                                      IN A GIVEN RANGE \n
	=======================================================================================\n"

case $?
in
0) spin=$(zenity --entry --text "ENTER THE SPIN YOU WANT (up/dn)")

case "$spin"
in
up)	vecmin=$(zenity --scale --text "pick your initial vector " --min-value=1 --max-value=30 --value=1 --step 2)
vecmax=$(zenity --scale --text "pick your final vector " --min-value=1 --max-value=30 --value=1 --step 2)
	(for vec in `seq ${vecmin} 1 ${vecmax}`;
	do
#	temp=`expr $vecmax - $vecmin`
#	temp=`expr 100 / $temp `
#	echo `expr $vec "*" $temp ` |zenity –progress –auto-close
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
echo "50";echo "# VECTOR NO . $vec GENERATED" ; sleep 0.2 
	cat >eigen_vec_spin${spin}${vec}.VASP_CHGCAR <vmdformat 
	awk -v vec=${vec} -f awk.awk <${eigenvecfile}${spin}.dat >>eigen_vec_spin${spin}${vec}.VASP_CHGCAR
	done;echo "100") |zenity  --pulsate --progress  --auto-close
	zenity --info --text " FILES GENERATED "  --timeout=2
	rm -f awk.awk info_temp ;;
dn) 	vecmin=$(zenity --scale --text "pick your initial vector " --min-value=1 --max-value=30 --value=1 --step 2)
vecmax=$(zenity --scale --text "pick your final vector " --min-value=1 --max-value=30 --value=1 --step 2)
	(for vec in `seq ${vecmin} 1 ${vecmax}`;
	do
#	temp=`expr $vecmax - $vecmin`
#	temp=`expr 100 / $temp `
#	echo `expr $vec "*" $temp ` |zenity –progress –auto-close
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
		}
		{i=col}{print $ i}
EOF
	cat >eigen_vec_spin${spin}${vec}.VASP_CHGCAR <vmdformat 
	echo "50";echo "# VECTOR NO . $vec GENERATED" ; sleep 0.2
	awk -v vec=${vec} -f awk.awk <${eigenvecfile}${spin}.dat >>eigen_vec_spin${spin}${vec}.VASP_CHGCAR
	echo "70" ;done;echo "100" )|zenity --progress --pulsate --auto-close
	zenity --info --text " FILES GENERATED " --timeout=2
	rm -f awk.awk info_temp ;;
*) inf=2 ;;
esac ;;

1) inf=2 ;;
esac
done
