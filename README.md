# Final Report and Project Software

## How to Run the Code
The below assumes that the user has `python` and `pip` already installed and accessible from the project directory.
Additionally, the commands below are for Windows users - Mac or Linux users should run the equivalent commands.

The easiest way to run the code is to first create and activate a virtual environment using
```
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```
Then install the required dependencies:

### Dependencies
Dependencies are listed in `requirements.txt`, and can be easily downloaded all at once using
```
pip install -r requirements.txt
```
The project requires `numba`, `numpy`, `matplotlib`, `ipykernel`, and `scipy`.

### Code Structure
Source code is structured as a Python module in `project/sandpiles`. 
To reproduce the experiment, the Jupyter notebook `project/sandpiles.ipynb` contains structured glue code and code to generate all figures and tables included in the report. It imports and calls functions defined in the `project/sandpiles` module.
Open the Jupyter notebook and select the virtual environment `.venv` just created as the kernel. You may need to run the code with `project` as your working directory.

