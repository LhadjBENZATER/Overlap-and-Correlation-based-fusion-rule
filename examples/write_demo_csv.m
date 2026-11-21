function write_demo_csv(path,header,values)
fid=fopen(path,'w'); assert(fid>=0,'Could not open CSV output');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',header);
for i=1:size(values,1)
 fprintf(fid,'%.10g',values(i,1));
 fprintf(fid,',%.10g',values(i,2:end)); fprintf(fid,'\n');
end
end
