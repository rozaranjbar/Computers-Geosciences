%% Generate_video

if do_SR

    video_name_SR = strcat('Video_SR_',datestr(now),'.avi'); 
    video_name_SR = strrep(video_name_SR,':','-');
    % create the video writer with 1 fps
    writerObj_SR = VideoWriter(video_name_SR);
    % set the seconds per image
    writerObj_SR.FrameRate = 10;  
    % open the video writer
    open(writerObj_SR);

    for t = 1:Tsim

        if t == 1

            real_time_SR_handle = figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        end 

        for i = 1:nx

            figure(real_time_SR_handle)
            subplot(3,5,i);
            plot(x_SR(i,1:t),'b','LineWidth',2)

            hold on
            box on
            grid on

            plot(upper_bound_SR(i,1:t),'b--','LineWidth',2)
            plot(lower_bound_SR(i,1:t),'b--','LineWidth',2)

            predicted_time = t:t+Np-1;
            plot(predicted_time+1,predicted_x_SR_selected(i,:,t),'r','LineWidth',2)
            plot(predicted_time,predicted_upper_bound_SR{selected_route(t)}(i,:),'r--','LineWidth',2)
            plot(predicted_time,predicted_lower_bound_SR{selected_route(t)}(i,:),'r--','LineWidth',2)
            title(strcat('x_{',num2str(i),'} in SR'))

            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';

            hold off

        end

        subplot(3,5,[11:15])
        plot([1:t],Visited_segment_SR(1:t),'bs--','MarkerFaceColor','b')
        hold on
        predicted_time = t:t+Np-1;
        plot(predicted_time,Route(selected_route(t),:),'rs--','MarkerFaceColor','r')
        hold off
        ylabel('time [s]')
        xlabel('visited segment')
        xticks(segments)
        title(strcat('Visited segment and expected route (t=',num2str(t),'s)'))

        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])

        % convert the image to a frame
        frame2write = getframe(gcf);
        % write the frames to the video
        writeVideo(writerObj_SR, frame2write);

    end 

end

if do_PD

    video_name_PD = strcat('Video_PD_',datestr(now),'.avi'); 
    video_name_PD = strrep(video_name_PD,':','-');
    % create the video writer with 1 fps
    writerObj_PD = VideoWriter(video_name_PD);
    % set the seconds per image
    writerObj_PD.FrameRate = 10; datestr(now) 
    % open the video writer
    open(writerObj_PD);

    for t = 1:Tsim

        if t == 1

            real_time_PD_handle = figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        end 

        for i = 1:nx

            figure(real_time_PD_handle)
            subplot(3,5,i);
            plot(x_PD(i,1:t),'b','LineWidth',2)

            hold on
            box on
            grid on

            plot(upper_bound_PD(i,1:t),'b--','LineWidth',2)
            plot(lower_bound_PD(i,1:t),'b--','LineWidth',2)

            predicted_time = t:t+Np-1;
            plot(predicted_time+1,predicted_x_PD(i,:,t),'r','LineWidth',2)
            plot(predicted_time,predicted_upper_bound_PD(i,:),'r--','LineWidth',2)
            plot(predicted_time,predicted_lower_bound_PD(i,:),'r--','LineWidth',2)
            title(strcat('x_{',num2str(i),'} in PD'))

            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';

            hold off

        end

        subplot(3,5,[11:15])
        plot([1:t],Visited_segment_PD(1:t),'bs--','MarkerFaceColor','b')
        hold on
        predicted_time = t:t+Np-1;
        Next_route_PD = Predifined_Route(t:t+Np-1); 
        plot(predicted_time,Next_route_PD,'rs--','MarkerFaceColor','r')
        hold off
        ylabel('time [s]')
        xlabel('visited segment')
        xticks(segments)
        title(strcat('Visited segment and expected route (t=',num2str(t),'s)'))

        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])

        % convert the image to a frame
        frame2write = getframe(gcf);
        % write the frames to the video
        writeVideo(writerObj_PD, frame2write);

    end 

end

% close the writer object
close(writerObj_SR);
close(writerObj_PD);
