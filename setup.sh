export DISABLE_AUTH=true
export ACCEPT_EULA=Y

export HTTP_PROXY=http://proxy01.pnl.gov:3128
export HTTPS_PROXY=http://proxy01.pnl.gov:3128
export NO_PROXY=*.pnl.gov,*.pnnl.gov,127.0.0.1,10.120.148.170,10.120.148.153

# To setup other dependencies
apt-get update -qq && apt-get install -y net-tools \
        python3.7 \
        python3-pip \
        curl \
        unixodbc

# To save to DB
curl https://packages.microsoft.com/ubuntu/20.10/prod/pool/main/m/msodbcsql17/msodbcsql17_17.7.2.1-1_amd64.deb --output mssql.deb
dpkg -i mssql.deb
curl https://packages.microsoft.com/config/ubuntu/20.04/prod.list | sudo tee /etc/apt/sources.list.d/msprod.list

Rscript -e "install.packages('argparse', repos='https://cran.microsoft.com/')"