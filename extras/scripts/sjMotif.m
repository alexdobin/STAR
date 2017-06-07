function [ ] = sjMotif( genomeDir, sjIn, sjOut )
%sjMotif add splice junction motif as the last column
%   genomeDir: STAR genome directory
%   sjIn:      splice junction loci file name (chr intron_strart intron_end)
%   sjOut:     output file name

%% load genome
fin1=fopen([genomeDir '/chrName.txt']);
chrName=textscan(fin1,'%s');
fclose(fin1);
chrName=chrName{1};

fin1=fopen([genomeDir '/chrStart.txt']);
chrStart=textscan(fin1,'%f');
fclose(fin1);
chrStart=chrStart{1};

fin1=fopen([genomeDir '/Genome']);
G=fread(fin1,inf,'*uint8');
fclose(fin1);

%% read sj file
fin1=fopen(sjIn);
sj=textscan(fin1,'%s %f %f %[^\n]');
fclose(fin1);

%% chr 
chrI=zeros(length(sj{1}),1);
for ii=1:length(chrName)
    chrI(strcmp(chrName{ii},sj{1}))=ii;
end

nT='ACGTN';
a=chrStart(chrI)+sj{2};
b=chrStart(chrI)+sj{3};
m=nT([G(a) G(a+1) G(b-1) G(b)]+1);

fout1=fopen(sjOut,'w');
for ii=1:length(chrI)
    fprintf(fout1,'%s\t%i\t%i\t%s\t%s\n',sj{1}{ii},sj{2}(ii),sj{3}(ii),sj{4}{ii},m(ii,:));
end
fclose(fout1);

end

