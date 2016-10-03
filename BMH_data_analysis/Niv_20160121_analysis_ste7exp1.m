function [] = Niv_20160121_analysis_ste7exp1()

%profile off
ipdir = 'C:\Users\Ben\Documents\GitHub\image_analysis\'
%adds image analysis directory to path based on directry of day's analysis
%code
path(ipdir,path)




%%
base_dir = 'C:\Users\Ben\Documents\Data\Niv\22_01_16Ste7exp1\'


% fname_save = 'processed_data_CD_1to3.mat';
% load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')
% 
% titles = {'C','D'}
% t = 1
% for jj = [1:2]
%     figure(1)
%     subplot(1,2,jj)
%     Ncells = length(all_tracks_vec{jj}.Exp)
%     nf_vec = []
%     nmi_vec = []
%     for kk = 1:Ncells
%         times = all_tracks_vec{jj}.Exp(kk).times
%         if ismember(t,times)
%             t_ind = find(times == t);
%             nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
%             nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
%         end
%     end
%     plot(nf_vec,nmi_vec,'x')
%     xlabel('nf vec')
%     ylabel('nmi vec')
%     title(titles{jj})
%     axis([0 16 0 11])
% end
% 
% all_tracks_vec;
% all_times_vec;



fname_save = 'processed_data_CD_allpos.mat';
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

titles = {'C','D'}
t = 1
for jj = [1:2]
    Ncells = length(all_tracks_vec{jj}.Exp);
    nf_vec = [];
    nmi_vec = [];
    for kk = 1:Ncells
        times = all_tracks_vec{jj}.Exp(kk).times;
        if ismember(t,times)
            t_ind = find(times == t);
            nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
            nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
        end
    end
    length(nf_vec)
    figure(2)
    subplot(1,2,jj)
    [X,bins] = hist3([nmi_vec',nf_vec'],[20,20]);
    imagesc(X)
    xlabel('nf vec')
    ylabel('nmi vec')
    title(titles{jj})
        
    
    figure(3)
    subplot(1,2,jj)
    plot(nf_vec,nmi_vec,'x')
    xlabel('nf vec')
    ylabel('nmi vec')
    title(titles{jj})
    axis([0 20 0 20])
    
    
    figure(4)
    subplot(1,2,jj)
    hist(nf_vec, 30)
    title('NF')
    xlim([0,25])
    
    figure(5)
    subplot(1,2,jj)
    hist(nmi_vec, 30)
    title('NMI')
    xlim([0,25])
end

% 
% %%Make Plot for all time points
% 
% titles = {'WT','IRE1 Del'}
% for t = [1:6]
%     for jj = [1:2]
%         Ncells = length(all_tracks_vec{jj}.Exp)
%         nf_vec = [];
%         nmi_vec = [];
%         for kk = 1:Ncells
%             times = all_tracks_vec{jj}.Exp(kk).times
%             if ismember(t,times)
%                 t_ind = find(times == t);
%                 nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
%                 nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
%             end
%         end
%         figure(6)
%         %2d Histogram (heatmap)
%         subplot(6,2, 2*(t-1) + jj)
%         [X,bins] = hist3([nmi_vec',nf_vec'],[20,20]);
%         imagesc(X)
%         xlabel('nf vec')
%         ylabel('nmi vec')
%         title([titles{jj}, ' time ',int2str(t-1)])
% 
% 
%         figure(7)
%         %2D Scatter
%         subplot(6,2,2*(t-1) + jj)
%         plot(nf_vec,nmi_vec,'x')
%         xlabel('nf vec')
%         ylabel('nmi vec')
%         title([titles{jj}, ' time ',int2str(t-1)])
%         axis([0 20 0 20])
% 
% 
%         figure(8)
%         %NF Histogram
%         subplot(6,2,2*(t-1) + jj)
%         hist(nf_vec, 30)
%         title(['NF ', titles{jj}, ' time ',int2str(t-1)])
%         xlim([0,25])
% 
%         figure(9)
%         %NMI Histogram
%         subplot(6,2,2*(t-1) + jj)
%         hist(nmi_vec, 30)
%         title(['NMI ', titles{jj}, ' time ',int2str(t-1)])
%         xlim([0,25])
%     end
% end
% 
% 
% 




return