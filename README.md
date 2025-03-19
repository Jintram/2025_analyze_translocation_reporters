

## Running the script

Script can be called as follows from the command line:

```bash
input_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Thrombine/"
output_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/Thrombine/"
auto_correct_bg=1

python analyze_transl_rep.py $input_folder $output_folder $auto_correct_bg nucleus 0 ERK 1 PKA 2
```

- `auto_correct_bg` should be 1 (=yes) or 0 (=no), indicating whether background should be corrected automatically.
- The keywords `nucleus 0 ERK 1 PKA 2` indicate that nuclear channel is 0, ERK measurement channel is 1, PKA measurement channel is 2. The keyword **'nucleus' is obligatory**, the rest can be changed as desired. At least one channel additional to the nuclear channel should be defined.
