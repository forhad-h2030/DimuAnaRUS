DIR_TOP=$(dirname $(readlink -f $BASH_SOURCE))
DIR_BUILD=$DIR_TOP/build
DIR_INST=$DIR_TOP/inst

#E1039_SHARE_VERSION=20240317
#E1039_CORE_VERSION=pr.155
#E1039_CORE_DIR=/seaquest/users/$USER/e1039/git/core-inst
FN_SETUP=/exp/seaquest/app/software/osg/software/e1039/this-e1039.sh
#FN_SETUP=/seaquest/users/mhossain/i_o_Manager/core-inst/this-e1039.sh
#FN_SETUP=/seaquest/users/mhossain/delete/core-inst/this-e1039.sh
if [ ! -e $FN_SETUP ] ; then # On grid
    FN_SETUP=/cvmfs/seaquest.opensciencegrid.org/seaquest/${FN_SETUP#/exp/seaquest/app/software/osg/}
fi
source $FN_SETUP
export   LD_LIBRARY_PATH=$DIR_INST/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$DIR_INST/include:$ROOT_INCLUDE_PATH

function cmake-this {
    if [ -e $DIR_BUILD ] ; then
	echo "Clean up the build directory..."
	\rm -rf $DIR_BUILD
    fi
    if [ -e $DIR_INST ] ; then
	echo "Clean up the install directory..."
	\rm -rf $DIR_INST
    fi
    mkdir -p $DIR_BUILD
    echo "Run cmake..."
    ( cd $DIR_BUILD && cmake -DCMAKE_INSTALL_PREFIX=$DIR_INST $DIR_TOP/src )
    local RET=$?
    if [ $RET -eq 0 ] ; then
	echo "OK.  Run 'make-this' at any time/directory."
    else
	echo "NG.  Need a fix."
    fi
    return $RET
}

function make-this {
    if [ ! -e $DIR_BUILD/Makefile ] ; then
	echo "The build directory is not ready.  Complete 'cmake-this' first."
	return 1
    fi
    ( cd $DIR_BUILD && make install )
    return $?
}

alias jobsub_q_mine='jobsub_q --group spinquest --user=$USER'
