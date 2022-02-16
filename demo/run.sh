echo -e "Generate Phantom cylinder!\n"

pushd ../python
python3 create_ball.py
popd


echo -e "\nCompiling MBIR Helical code!\n"
pushd ../src
make
popd 
