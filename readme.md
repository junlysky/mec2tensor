# 说明

## mec2info改编自focmec中的dsretc.f，使用了fortran90标准，实现了批量处理的功能。
## 程序功能为由震源机制解一个断层面的strike,dip,rake计算地震矩张量的全部信息，震源模型为双力偶模型。
## 输入文件名词必须为：mec.d
## mec.d的格式必须为：lon lat strike dip rake
## 输出文件为：mecinfo.dat  
