import subprocess
p = subprocess.Popen("G:\Programmes\GPT\\bin\gpt.exe", cwd=r"G:\Programmes\GPT\bin")
stdout, stderr = p.communicate()