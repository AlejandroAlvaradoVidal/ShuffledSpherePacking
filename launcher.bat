@echo off
:: 1. Activar tu entorno de Conda
call C:\Users\ale_m\miniconda3\Scripts\activate.bat MSFA

:: 2. Ejecutar el código y guardar todo lo que saldría en la consola en un archivo .txt
python PythonCode/gen_cssp_v4_nf31.py > log_progreso.txt 2>&1