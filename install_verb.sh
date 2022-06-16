# create ~/.local/bin if it does not exist
local_bin="$HOME/.local/bin"
if [ ! -d $local_bin ]
then
	echo "creating folder $local_bin"
	mkdir $local_bin
fi

# create ./build folder if it does not exist
if [ ! -d ./build ]
then
    mkdir ./build
fi

# run cmake and make to compile VERB4D_Solver and install it to ~/.local/bin
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$HOME/.local"
make -j2
make install
cd ..

# if ~/.local/bin is not part of the PATH variable append the definition to ~/.bashrc
case ":$PATH:" in
    *:$local_bin*)
    ;;
    *)
	echo "adding $local_bin to PATH in .bashrc"
	export PATH=$local_bin:$PATH
	echo "export PATH=$local_bin:"'$PATH' >> ~/.bashrc
    ;;
esac
