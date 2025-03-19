

## Running the script

Script can be called as follows from the command line:

```bash
input_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Thrombine/"
output_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/"
auto_correct_background=1

python analyze_transl_rep.py $input_folder $output_folder $auto_correct_background nucleus 0 ERK 1 PKA 2
```

The keywords `nucleus 0 ERK 1 PKA 2` indicate that nuclear channel is 0, ERK measurement channel is 1, PKA measurement channel is 2.