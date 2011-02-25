#!/bin/bash

INSTALL_DIR=
MPI_PATH=
CONF_FILE=

die() {
	echo $@
	exit 1
}

for x in $@; do
	case $x in
		--with-mpi=*)
			MPI_PATH=${x#--with-mpi=}
			;;
		--prefix=*)
			INSTALL_DIR=${x#--prefix=}
			;;
		--conf=*)
			CONF_FILE=${x#--conf=}
			;;
	esac
done

[[ -z $OMD_HOME ]] || INSTALL_DIR=$OMD_HOME
[[ -z $INSTALL_DIR ]] && die 'install target directory needed as $1'
[[ -z $CONF_FILE ]] && CONF_FILE=$HOME/.omdconf-$(uname -m)

echo using target directory: $INSTALL_DIR

if [[ -d $INSTALL_DIR ]]; then
	echo removing existing target directory..
	echo replaced by the new version.
	rm -rf $INSTALL_DIR
fi

mkdir -p \
	$INSTALL_DIR/include \
	$INSTALL_DIR/lib \
	$INSTALL_DIR/tables

cp -r include/* $INSTALL_DIR/include
cp -r lib/* $INSTALL_DIR/lib
cp -r tables/* $INSTALL_DIR/tables

if [[ -f $CONF_FILE ]]; then
	echo the configuration file exist. keeping the file...
else
cat << endl > $CONF_FILE
export OMD_HOME=$INSTALL_DIR
export OMD_TABLE=\$OMD_HOME/tables
export OMD_CLASS=\$OMD_HOME/include/class
export OMD_LIBS="-I$OMD_HOME/include -L$OMD_HOME/lib -lomd-$(uname -m) -lomd_object-$(uname -m) -lm"
endl

[[ ! -z $MPI_PATH ]] && echo "export OMD_MPI=$MPI_PATH" >> $CONF_FILE

fi

exit 0
