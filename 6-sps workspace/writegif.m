function res=writegif(name,frames,dt)
nframe=length(frames);
for i=1:nframe
    [image,map]=frame2im(frames(i));
    [im,map2]=rgb2ind(image,128);
    if i==1
        imwrite(im,map2,name,'gif','writeMode','overwrite','delaytime',dt,'loopcount',inf);
    else
        imwrite(im,map2,name,'writeMode','append','delaytime',dt); %,'loopcount',inf);
    end
end
end