MODEL =  1;

%预留切割机切口宽度：0.1mm

%原构型
if MODEL == 0
    sketch_scr = Z_cal_scr( P(BEST_NUM(2),:,BEST_NUM(1)) ) ;
    SHAPE = '-Z';
end

%加宽构型
if MODEL == 1
    tar.range(1,:) = [ 1, 1.5];
    sketch_scr = Z_cal_scr_W( P(BEST_NUM(2),:,BEST_NUM(1)) ) ;
    SHAPE = '-H';
end


for k = 1 : sketch_scr.num_line
    
    line1 = sketch_scr.line1{k};
    line2 = flip( sketch_scr.line2{k} );
    line3 = sketch_scr.line3{k};
    line4 = sketch_scr.line4{k};
    
    name_scr = ['Cylinder-sketch-',ORDER,SHAPE,num2str(k),'.scr'];
    
    fid=fopen(name_scr,'w');
    fprintf(fid,'line\n');
    for i = 1:length(line1)
        fprintf(fid,'%g,%g\n',line1(i,1),line1(i,2));
    end

    for i = 1:length(line4)
        fprintf(fid,'%g,%g\n',line4(i,1),line4(i,2));
    end

    for i = 1:length(line2)
        fprintf(fid,'%g,%g\n',line2(i,1),line2(i,2));
    end

    for i = 1:length(line4)
        fprintf(fid,'%g,%g\n',line3(i,1),line3(i,2));
    end

    fclose(fid);
    
end




