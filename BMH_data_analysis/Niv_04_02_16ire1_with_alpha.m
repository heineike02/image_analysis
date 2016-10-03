function [] = Niv_04_02_16ire1_with_alpha()

%profile off

%%
base_dir = 'C:\Users\Ben\Documents\Data\Niv\04_02_16ire1_with_alpha\'


fname_save = 'processed_data'
load([base_dir,fname_save],'all_times_vec','all_tracks_vec','posvec')

%This is the well order for the dataset: 
%'A3','B3','C3','D3','E3','F3','G3','H3','G2','H2'


% titles = {'WT 5uM alpha','ire1 del 5uM alpha','C3','D3','E3','F3','G3','H3','G2','H2'}
% t = 1
% for jj = [1:length(titles)]
%     figure(1)
%     Ncells = length(all_tracks_vec{jj}.Exp);
%     nf_vec = [];
%     nmi_vec = [];
%     for kk = 1:Ncells
%         times = all_tracks_vec{jj}.Exp(kk).times;
%         if ismember(t,times)
%             t_ind = find(times == t);
%             nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
%             nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
%         end
%     end
%     
%     figure(1)
%     %2d Histogram (heatmap)
%     subplot(3,4,jj)
%     length(nf_vec)
%     [X,bins] = hist3([nmi_vec',nf_vec'],[20,20]);
%     imagesc(X)
%     xlabel('nf vec')
%     ylabel('mi vec')
%     title(titles{jj})
%         
%     
%     figure(2)
%     %Scatter plot
%     subplot(3,4,jj)
%     plot(nf_vec,nmi_vec,'x')
%     xlabel('nf vec')
%     ylabel('mi vec')
%     title(titles{jj})
%     axis([0 20 0 20])
%     
%     
%     figure(3)
%     %NF Histogram
%     subplot(3,4,jj)
%     hist(nf_vec, 30)
%     title(['NF ', titles{jj}])
%     xlim([0,25])
%     
%     figure(4)
%     %Median Intensity Histogram
%     subplot(3,4,jj)
%     hist(nmi_vec, 30)
%     title(['MI ', titles{jj}])
%     xlim([0,25])
% end

% %%Make Plot for all time points for the first 6 wells
% titles = {'A3','B3','C3','D3','E3','F3'}
% ntimes = 7
% for t = [1:ntimes]
%     for jj = [1:length(titles)]
%         Ncells = length(all_tracks_vec{jj}.Exp);
%         nf_vec = [];
%         nmi_vec = [];
%         for kk = 1:Ncells
%             times = all_tracks_vec{jj}.Exp(kk).times;
%             if ismember(t,times)
%                 t_ind = find(times == t);
%                 nf_vec = [nf_vec, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_ind)];
%                 nmi_vec = [nmi_vec, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_ind)];
%             end
%         end
%         
%         time_val = num2str(all_times_vec{jj}.Exp(t));
%         
%         figure(6)
%         %2d Histogram (heatmap)
%         subplot(length(titles),ntimes, (jj-1)*ntimes + t)
%         [X,bins] = hist3([nmi_vec',nf_vec'],[20,20]);
%         imagesc(X)
%         xlabel('nf vec')
%         ylabel('nmi vec')
%         title([titles{jj}, ' time = ',time_val])
% 
% 
%         figure(7)
%         %2D Scatter
%         subplot(length(titles),ntimes, (jj-1)*ntimes + t)
%         plot(nf_vec,nmi_vec,'x')
%         xlabel('nf vec')
%         ylabel('nmi vec')
%         title([titles{jj}, ' time = ',time_val])
%         axis([0 20 0 20])
% 
% 
%         figure(8)
%         %NF Histogram
%         subplot(length(titles),ntimes, (jj-1)*ntimes + t)
%         hist(nf_vec, 30)
%         title(['NF ', titles{jj}, ' time = ',time_val])
%         xlim([0,25])
% 
%         figure(9)
%         %NMI Histogram
%         subplot(length(titles),ntimes, (jj-1)*ntimes + t)
%         hist(nmi_vec, 30)
%         title(['NMI ', titles{jj}, ' time = ',time_val])
%         xlim([0,25])
%     end
% end

%%Make delta nf plots for each condition
titles = {'WT 5uM alpha','ire1 del 5uM alpha', 'WT 10uM alpha','ire1 del 10uM alpha', 'WT','ire1 del'}

start_and_stop_time = [1,7];

for jj = [1:length(titles)]
    Ncells = length(all_tracks_vec{jj}.Exp);
    nf_vec_start = [];
    nmi_vec_start = [];
    nf_vec_stop = [];
    nmi_vec_stop = [];
    
    for kk = 1:Ncells
        times = all_tracks_vec{jj}.Exp(kk).times;
        if sum(ismember(start_and_stop_time,times)) == 2
            t_start_ind = find(times == start_and_stop_time(1));
            t_stop_ind = find(times == start_and_stop_time(2));
            nf_vec_start = [nf_vec_start, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_start_ind)];
            nf_vec_stop = [nf_vec_stop, all_tracks_vec{jj}.Exp(kk).nf.RFP(t_stop_ind)];
            nmi_vec_start = [nmi_vec_start, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_start_ind)];
            nmi_vec_stop = [nmi_vec_stop, all_tracks_vec{jj}.Exp(kk).nmi.RFP(t_stop_ind)];
        end
    end

    

%     figure(6)
%     %2d Histogram (heatmap)
%     subplot(length(titles),ntimes, (jj-1)*ntimes + t)
%     [X,bins] = hist3([nmi_vec',nf_vec'],[20,20]);
%     imagesc(X)
%     xlabel('nf vec')
%     ylabel('nmi vec')
%     title([titles{jj}, ' time = ',time_val])

    length(nf_vec_start)
    figure(10)
    %2D Scatter
    subplot(3,2,jj)
    plot(nf_vec_start, nf_vec_start-nf_vec_stop ,'x')
    xlabel('Initial Nuclear Fluorescence Score')
    ylabel('Initial - Final Nuclear Fluorescence Score')
    title(titles{jj})
    axis([0 35 -20 30])


    figure(11)
    %NF Histogram
    subplot(3,2,jj)
    hist(nf_vec_start - nf_vec_stop, 30)
    xlabel('Initial - Final Nuclear Fluorescence Score')
    title(['NF ', titles{jj}])
    xlim([-20,30])
% 
%     figure(9)
%     %NMI Histogram
%     subplot(length(titles),ntimes, (jj-1)*ntimes + t)
%     hist(nmi_vec, 30)
%     title(['NMI ', titles{jj}, ' time = ',time_val])
%     xlim([0,25])
end







return