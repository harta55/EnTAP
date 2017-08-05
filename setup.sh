#TODO update setup script
echo "-----------------------"
echo "Setting up DIAMOND"
echo "-----------------------"
cd libs/diamond-0.8.31/ || \
	{ echo "libs/diamond-0.8.31 path not found"; exit 1; }
mkdir -p bin
cd bin
(cmake .. &&  make) || \
	{ echo "Error in compiling DIAMOND"; exit 1; }
echo "DIAMOND setup complete!"
echo "-----------------------"
echo "Setting up RSEM"
echo "-----------------------"
cd ../../RSEM-1.3.0 || \
	{ echo "libs/RSEM-1.3.0 path not found"; exit 1; }
(make && make  ebseq) || \
	{ echo "Error in compiling RSEM"; exit 1; }
echo "Success! DIAMOND and RSEM compiled, continue to EnTAP"
