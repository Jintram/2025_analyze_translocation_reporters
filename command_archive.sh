




# If you store how you executed your script in a text file, 
# you can easily find back what you did later.


###

input_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Thrombine/"
output_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/Thrombine/"
auto_correct_bg=1

python analyze_transl_rep.py $input_folder $output_folder $auto_correct_bg nucleus 0 ERK 1 PKA 2



###
input_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_DATA_julian/Forskolin/"
output_folder="/Users/m.wehrens/Data_UVA/2024_10_Sebastian-KTR/202503_OUTPUT-testmw/Forskolin/"
auto_correct_bg=1

python analyze_transl_rep.py $input_folder $output_folder $auto_correct_bg nucleus 0 ERK 1 PKA 2