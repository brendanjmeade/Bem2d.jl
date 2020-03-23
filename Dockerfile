FROM julia:latest

# a fortran compiler is necessary for installing okada_wrapper
RUN apt-get update && apt-get install -y gfortran

COPY . /app
WORKDIR /app

# This will install all the dependencies specified in Project.toml
RUN julia -e 'import Pkg;\
    Pkg.add(Pkg.PackageSpec(path="."));\
    Pkg.add("IJulia");'

# equivalent of `pip install okada_wrapper`
RUN julia -e 'import Pkg;\
    Pkg.activate(".");\
    using PyCall;\
    pyimport("pip").main(["install", "okada_wrapper"])\
    '

# precompile all the packages we installed so that it's fast to get up and
# running once we start a container
RUN julia -e 'import Pkg; Pkg.activate("."); Pkg.precompile();'

# finally, use conda to install jupyterlab
RUN /root/.julia/conda/3/bin/conda install -y jupyterlab

CMD ["bash"]
