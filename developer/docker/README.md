# More information about NMRlipids Docker images 

Right now there are:

1. `nmrlipids/gromacs`: This is the base gromacs image, only contains gromacs.

Then we have one image extending this gromacs base:

3. `nmrlipids/core`: features everything developers need to work on the project 


# Installation using Docker

We provide a Docker-based development environment that allows for easy testing of new code and features. This method is recommended for development and testing purposes.

## 1. Install Docker

Before using the Docker-based development environment, you'll need to have Docker installed on your system:

- **Linux**: Follow the official Docker installation guide for Linux: [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/)
- **macOS**: Download and install Docker Desktop for Mac from [https://www.docker.com/products/docker-desktop/](https://docs.docker.com/desktop/setup/install/mac-install/)
- **Windows**: Download and install Docker Desktop for Windows from [https://www.docker.com/products/docker-desktop/](https://docs.docker.com/desktop/setup/install/windows-install/) (requires Windows Subsystem for Linux)

## 2. Using the Docker Development Environment

### Setting Up the Development Environment
Go to the directory where you have the fairmd-lipids repository. Then,

1. Download the latest fairmd-lipids core image:
   ```
   docker pull fairmd-lipids/core:latest
   ```

2. Alternatively the Docker image can be built locally:
   ```bash
   docker build -t NAME_OF_THE_DOCKER_IMAGE .
   ```


2. Initialize the submodule data:
   ```bash
   git submodule update --init
   ```

### Testing Code

1. Start the container with your code mounted:
   ```bash
   docker run -it -v $(pwd):/github/workspace NAME_OF_THE_DOCKER_IMAGE
   ```
   or
   if you pulled the latest image,
   ```bash
   docker run -it -v $(pwd):/github/workspace fairmd-lipids/core:latest
   ```

2. Inside the container:
   ```bash
   # Install base requirements and the fairmd-lipids package in editable mode
   pip install -e .
   
   # Run tests
   ./runtests.sh
   ./runtests-regression.sh
   ```

   You can specify additional arguments to the `pytest` via
   `tox -e tests`

   For full list of accepted flags, please look at [pytest documentation](https://docs.pytest.org/en/6.2.x/usage.html).

   Do not forget to add this new package to the requirements.txt! 

3. When done testing:
   ```bash
   # Exit container
   exit
   ```

### Managing Different Test Environments

By default the core image features the dev requirements from this repository but you can easily create different testing environments with different dependencies:

1. Start container again:
   ```bash
   docker run -it -v $(pwd):/workspace NAME_OF_THE_DOCKER_IMAGE
   ```

2. Inside container, create a new virtual environment using Conda as previously described. (Must update this)

### Tips for Docker Development

- Your local code changes are automatically reflected in the container since we use a volume mount
- You can run multiple test environments simultaneously by creating different virtual environments
- To see all running containers: `docker ps`
- To stop a container: `docker stop <container_id>`
- To remove all resources related to docker including images, container, build-resources: ` docker system prune -a`
