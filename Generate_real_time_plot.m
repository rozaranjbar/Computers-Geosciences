%% Generate_real_time_plot


        if video

            
    
        end


    if do_SR
        
        if video

            % open the video writer
            open(writerObj_SR);
    
        end

        if t == 1

            real_time_SR_handle = figure;
            
            if ~video
            
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0 0.5 1]);
            
            else
                
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                
            end

        end 
        
        if segments_used == 2

            for i = 1:nx

                figure(real_time_SR_handle)
                subplot(6,6,i);
                
                plot(x_controller_SR(i,:),'b','LineWidth',2)
                
                hold on
                
                plot(x_real_SR(i,:),'b','LineWidth',2)
                box on
                grid on

                plot(upper_bound_SR(i,:),'b--','LineWidth',2)
                plot(lower_bound_SR(i,:),'b--','LineWidth',2)

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

            subplot(6,6,[31:36])
            plot(Visited_segment_SR(1:t),[1:t],'bs--','MarkerFaceColor','b')
            hold on
            predicted_time = t:t+Np-1;
            plot(Route(selected_route(t),:),predicted_time,'rs--','MarkerFaceColor','r')
            hold off
            ylabel('time [s]')
            xlabel('visited segment')
            xticks(segments)
            title(strcat('Visited segment and expected route (t=',num2str(t),'s)'))

            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([1 nx 0 Tsim])
        
        elseif segments_used == 1
            
            for i = 1:ns

                figure(real_time_SR_handle)
                subplot(3,4,i);
                
                plot(x_controller_SR(segments(i),:),'b','LineWidth',2)
                
                hold on
                
                plot(x_real_SR(segments(i),:),'g','LineWidth',2)

                box on
                grid on

                plot(upper_bound_SR(segments(i),:),'b--','LineWidth',2)
                plot(lower_bound_SR(segments(i),:),'b--','LineWidth',2)

                predicted_time = t:t+Np-1;
                plot(predicted_time+1,predicted_x_SR_selected(segments(i),:,t),'r','LineWidth',2)
                plot(predicted_time,predicted_upper_bound_SR{selected_route(t)}(segments(i),:),'r--','LineWidth',2)
                plot(predicted_time,predicted_lower_bound_SR{selected_route(t)}(segments(i),:),'r--','LineWidth',2)
                title(strcat('x_{',num2str(segments(i)),'} in SR'))

                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';

                hold off

            end

            subplot(3,4,[9:12])
            plot(Visited_segment_SR(1:t),[1:t],'bs--','MarkerFaceColor','b')
            hold on
            predicted_time = t:t+Np-1;
            plot(Route(selected_route(t),:),predicted_time,'rs--','MarkerFaceColor','r')
            hold off
            ylabel('time [s]')
            xlabel('visited segment')
            xticks(segments)
            title(strcat('Visited segment and expected route (t=',num2str(t),'s)'))

            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([1 nx 0 Tsim])
            
        end
        
        if video

            % convert the image to a frame
            frame2write = getframe(gcf);
            % write the frames to the video
            writeVideo(writerObj_SR, frame2write);
    
        end

    end 

    if do_PD
        
        if video

            % open the video writer
            open(writerObj_PD);
    
        end

        if t == 1

            real_time_PD_handle = figure;
            
            if ~video
            
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 1]);
            
            else
                
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                
            end

        end 
        
        if segments_used == 2

            for i = 1:nx

                figure(real_time_PD_handle)
                subplot(6,6,i);
                
                plot(x_controller_PD(i,:),'b','LineWidth',2)
                
                hold on
                
                plot(x_real_PD(i,:),'g','LineWidth',2)

                box on
                grid on

                plot(upper_bound_PD(i,:),'b--','LineWidth',2)
                plot(lower_bound_PD(i,:),'b--','LineWidth',2)

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

            subplot(6,6,[31:36])
            plot(Visited_segment_PD(1:t),[1:t],'bs--','MarkerFaceColor','b')
            hold on
            predicted_time = t:t+Np-1;
            Next_route_PD = Predefined_Route(t:t+Np-1); 
            plot(Next_route_PD,predicted_time,'rs--','MarkerFaceColor','r')
            hold off
            ylabel('time [s]')
            xlabel('visited segment')
            xticks(segments)
            title(strcat('Visited segment and expected route (t=',num2str(t),'s)'))

            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([1 nx 0 Tsim])
        
        elseif segments_used == 1
            
            for i = 1:ns

                figure(real_time_PD_handle)
                subplot(3,4,i);
                
                plot(x_controller_PD(segments(i),:),'b','LineWidth',2)
                
                hold on
                
                plot(x_real_PD(segments(i),:),'g','LineWidth',2)
       
                box on
                grid on

                plot(upper_bound_PD(segments(i),:),'b--','LineWidth',2)
                plot(lower_bound_PD(segments(i),:),'b--','LineWidth',2)

                predicted_time = t:t+Np-1;
                plot(predicted_time+1,predicted_x_PD(segments(i),:,t),'r','LineWidth',2)
                plot(predicted_time,predicted_upper_bound_PD(segments(i),:),'r--','LineWidth',2)
                plot(predicted_time,predicted_lower_bound_PD(segments(i),:),'r--','LineWidth',2)
                title(strcat('x_{',num2str(segments(i)),'} in PD'))

                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';

                hold off

            end

            subplot(3,4,[9:12])
            plot(Visited_segment_PD(1:t),[1:t],'bs--','MarkerFaceColor','b')
            hold on
            predicted_time = t:t+Np-1;
            Next_route_PD = Predefined_Route(t:t+Np-1); 
            plot(Next_route_PD,predicted_time,'rs--','MarkerFaceColor','r')
            hold off
            ylabel('time [s]')
            xlabel('visited segment')
            xticks(segments)
            title(strcat('Visited segment and expected route (t=',num2str(t),'s)'))

            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([1 nx 0 Tsim])
            
        end
        
        if video

            % convert the image to a frame
            frame2write = getframe(gcf);
            % write the frames to the video
            writeVideo(writerObj_PD, frame2write);
    
        end

    end 