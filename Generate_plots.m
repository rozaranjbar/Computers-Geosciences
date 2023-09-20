%% Generate_plots

%% Legends

legends_x = ["x_{1}","x_{2}","x_{3}","x_{4}","x_{5}","x_{6}","x_{7}","x_{8}","x_{9}","x_{10}"];
legends_water_levels = ["x_{1}","x_{3}","x_{7}","x_{9}","x_{11}","x_{21}","x_{25}"];
legends_u = ["u_{1}","u_{2}","u_{3}","u_{4}","u_{5}","u_{6}","u_{7}","u_{8}"];

%%  

if individual_plots

    if do_SR

        figure

        subplot(2,1,1)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(x_real_SR(i,:),'LineWidth',2)

        end
        title('x in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_x)

        subplot(2,1,2)
        hold on
        box on
        grid on

        for i = 1:size(water_levels_SR,1)

            plot(water_levels_SR(i,:),'LineWidth',2)

        end

        title('water levels in SR')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_water_levels)

    end

end

%%

if individual_plots

    if do_PD

        figure

        subplot(2,1,1)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(x_real_PD(i,:),'LineWidth',2)

        end
        title('x in predefined route')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_x)

        subplot(2,1,2)
        hold on
        box on
        grid on

        for i = 1:size(water_levels_PD,1)

            plot(water_levels_PD(i,:),'LineWidth',2)

        end

        title('water levels in PD')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_water_levels)

    end

end

%%

if individual_plots

    if do_MPC

        figure
        subplot(2,1,1)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(x_real_MPC(i,:),'LineWidth',2)

        end
        title('x in classical MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_x)

        subplot(2,1,2)
        hold on
        box on
        grid on

        for i = 1:size(water_levels_MPC,1)

            plot(water_levels_MPC(i,:),'LineWidth',2)

        end

        title('water levels in classical MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_water_levels)

    end

end

%%

if individual_plots

    if do_SR

        figure
        hold on
        box on
        grid on

        for i = 1:nu

            plot(applied_u_SR(i,:),'LineWidth',2)

        end
        title('applied u in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_u)

    end

end

%%

if individual_plots

    if do_PD

        figure
        hold on
        box on
        grid on

        for i = 1:nu

            plot(applied_u_PD(i,:),'LineWidth',2)

        end
        title('applied u in predefined route')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_u)

    end

end

%%

if individual_plots

    if do_MPC

        figure
        hold on
        box on
        grid on

        for i = 1:nu

            plot(applied_u_MPC(i,:),'LineWidth',2)

        end
        title('applied u in classical MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_u)

    end

end


%%

if ~individual_plots

    if do_SR && do_PD && do_MPC

        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        subplot(3,3,1)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(x_real_SR(i,:),'LineWidth',2)
           
        end
        title('x in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_x)

        subplot(3,3,2)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(x_real_PD(i,:),'LineWidth',2)

        end
        title('x in predefined route')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_x)

        subplot(3,3,3)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(x_real_MPC(i,:),'LineWidth',2)

        end
        title('x in classical MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_x)

        subplot(3,3,7)
        hold on
        box on
        grid on

        for i = 1:nu

            plot(applied_u_SR(i,:),'LineWidth',2)

        end
        title('applied u in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_u)

        subplot(3,3,8)
        hold on
        box on
        grid on

        for i = 1:nu

            plot(applied_u_PD(i,:),'LineWidth',2)

        end
        title('applied u in predefined route')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_u)

        subplot(3,3,9)
        hold on
        box on
        grid on

        for i = 1:nu

            plot(applied_u_MPC(i,:),'LineWidth',2)

        end
        title('applied u in classical MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_u)

        subplot(3,3,4)
        hold on
        box on
        grid on

        for i = 1:size(water_levels_SR,1)

            plot(water_levels_SR(i,:),'LineWidth',2)

        end

        title('water levels in SR')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_water_levels)

        subplot(3,3,5)
        hold on
        box on
        grid on

        for i = 1:size(water_levels_PD,1)

            plot(water_levels_PD(i,:),'LineWidth',2)

        end

        title('water levels in PD')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_water_levels)

        subplot(3,3,6)
        hold on
        box on
        grid on

        for i = 1:size(water_levels_MPC,1)

            plot(water_levels_MPC(i,:),'LineWidth',2)

        end

        title('water levels in classical MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        legend(legends_water_levels)  

    end

end

%%

if individual_plots

    if do_SR

        figure
        hold on
        box on
        grid on

        for i = 1:nx

            plot(sigma_xi_SR(i,:),'LineWidth',2)

        end

        yyaxis right

        for i = 1:nx

            plot(water_level_mean_disturbances(i,:),'LineWidth',2,'LineStyle','--')

        end

        title('\sigma_{x} in proposed algorithm')
        legend('\sigma_{x1}','\sigma_{x2}','\sigma_{x3}','\sigma_{x4}','\sigma_{x5}','\sigma_{x6}','\sigma_{x7}','\sigma_{x8}','\sigma_{x9}','\sigma_{x10}',...
            '\mu_{x1}','\mu_{x2}','\mu_{x3}','\mu_{x4}','\mu_{x5}','\mu_{x6}','\mu_{x7}','\mu_{x8}','\mu_{x9}','\mu_{x10}')
        lgd.Layout.Tile = 'east';
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';

    end

end

%%

if individual_plots

    if do_PD

        figure
        hold on
        box on
        grid on

        for i = 1:nx

            plot(sigma_xi_PD(i,:),'LineWidth',2)

        end

        yyaxis right

        for i = 1:nx

            plot(water_level_mean_disturbances(i,:),'LineWidth',2,'LineStyle','--')

        end

        title('\sigma_{x} in predefined route')
        legend('\sigma_{x1}','\sigma_{x2}','\sigma_{x3}','\sigma_{x4}','\sigma_{x5}','\sigma_{x6}','\sigma_{x7}','\sigma_{x8}','\sigma_{x9}','\sigma_{x10}',...
            '\mu_{x1}','\mu_{x2}','\mu_{x3}','\mu_{x4}','\mu_{x5}','\mu_{x6}','\mu_{x7}','\mu_{x8}','\mu_{x9}','\mu_{x10}')
        lgd.Layout.Tile = 'east';
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
    end

end

%% 

if ~individual_plots

    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    if do_SR

        subplot(2,1,1)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(sigma_xi_SR(i,:),'LineWidth',2)

        end

        yyaxis right

        for i = 1:nx

            plot(water_level_mean_disturbances(i,:),'LineWidth',2,'LineStyle','--')

        end

%             title('\sigma_{x} in proposed algorithm')
%             legend('\sigma_{x1}','\sigma_{x2}','\sigma_{x3}','\sigma_{x4}','\sigma_{x5}','\sigma_{x6}','\sigma_{x7}','\sigma_{x8}','\sigma_{x9}','\sigma_{x10}',...
%                 '\mu_{x1}','\mu_{x2}','\mu_{x3}','\mu_{x4}','\mu_{x5}','\mu_{x6}','\mu_{x7}','\mu_{x8}','\mu_{x9}','\mu_{x10}')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';

    end

    if do_PD

        subplot(2,1,2)
        hold on
        box on
        grid on

        for i = 1:nx

            plot(sigma_xi_PD(i,:),'LineWidth',2)

        end

        yyaxis right

        for i = 1:nx

            plot(water_level_mean_disturbances(i,:),'LineWidth',2,'LineStyle','--')

        end

        title('\sigma_{x} in predefined route')
        legend('\sigma_{x1}','\sigma_{x2}','\sigma_{x3}','\sigma_{x4}','\sigma_{x5}','\sigma_{x6}','\sigma_{x7}','\sigma_{x8}','\sigma_{x9}','\sigma_{x10}',...
            '\mu_{x1}','\mu_{x2}','\mu_{x3}','\mu_{x4}','\mu_{x5}','\mu_{x6}','\mu_{x7}','\mu_{x8}','\mu_{x9}','\mu_{x10}')
        lgd.Layout.Tile = 'east';
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';

    end

end

%%

if do_SR && do_PD && do_MPC

    figure

    hold on
    box on
    grid on

    for i = 1:3

        plot(comparison(i,:),'LineWidth',2)

    end
    legend('Proposed Algorithm','Predifined Route','Classical MPC')
    title('Comparison stage cost')
    ax = gca;
    ax.FontSize = 14;
    ax.FontWeight = 'bold';

end

%%

if individual_plots

    if do_SR

        figure
        hold on
        box on
        grid on

        plot(Visited_segment_SR(1:last_t),[1:last_t],'rs--','MarkerFaceColor','r')
%         plot([1:last_t],Visited_segment_SR(1:last_t),'rs--','MarkerFaceColor','r')

        ylabel('time [step]')
        xlabel('visited segments')
        xticks(segments)
        title('Visited segments in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])

    end

end

%% 

if ~individual_plots

    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    if do_SR

        subplot(2,2,1)
        hold on
        box on
        grid on

        plot(Visited_segment_SR(1:last_t),[1:last_t],'rs--','MarkerFaceColor','r')
%         plot([1:last_t],Visited_segment_SR(1:last_t),'rs--','MarkerFaceColor','r')

        ylabel('time [step]')
        xlabel('visited segments')
        xticks(segments)
        title('Visited segments in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])
        
        subplot(2,2,2)
        hold on
        box on
        grid on

        plot([1:last_t],robot_battery_SR(1:last_t),'r-','MarkerFaceColor','r','LineWidth',2)

        xlabel('time [step]')
        ylabel('Robot battery (%)')
        title('Robot battery in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])


    end

    if do_PD

        subplot(2,2,3)
        hold on
        box on
        grid on

        plot(Predefined_Route(1:last_t),[1:last_t],'rs--','MarkerFaceColor','r')
%         plot([1:last_t],Predefined_Route(1:last_t),'rs--','MarkerFaceColor','r')

        ylabel('time [step]')
        xlabel('visited segments')
        xticks(segments)
        title('Visited segments in PR-MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])
        
        subplot(2,2,4)
        hold on
        box on
        grid on

        plot([1:last_t],robot_battery_PD(1:last_t),'r-','MarkerFaceColor','r','LineWidth',2)

        xlabel('time [step]')
        ylabel('Robot battery (%)')
        title('Robot battery in PR-MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])

    end

end

%% 

if individual_plots

    if do_SR

        figure
        hold on
        box on
        grid on

        plot([1:last_t],robot_battery_SR(1:last_t),'r-','MarkerFaceColor','r')

        xlabel('time [step]')
        ylabel('Robot battery (%)')
        title('Robot battery in proposed algorithm')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])

    end

    if do_PD

        figure
        hold on
        box on
        grid on

        plot([1:last_t],robot_battery_PD(1:last_t),'r-','MarkerFaceColor','r')

        xlabel('time [step]')
        ylabel('Robot battery (%)')
        title('Robot battery in PR-MPC')
        ax = gca;
        ax.FontSize = 14;
        ax.FontWeight = 'bold';
        axis([1 nx 0 Tsim])

    end

end

%%

if ~individual_plots
    
    if nx <= 10

        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        index_subplot = 1;

        if do_SR

            for i = 1:nx

                variable_to_plot = i;

                for t = 1:Tsim

                    upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_SR(variable_to_plot,t),sat_max_constraint);
                    lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_SR(variable_to_plot,t),sat_min_constraint);

                end

                if index_subplot <= 5

                    subplot(4,5,index_subplot)

                else

                   subplot(4,5,index_subplot+5) 

                end

                hold on
                box on
                grid on

                plot(x_SR(variable_to_plot,:),'b','LineWidth',2)
                plot(upper_bound,'b--','LineWidth',2)
                plot(lower_bound,'b--','LineWidth',2)

                title(strcat('x_{',num2str(i),'} in proposed algorithm'))
                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';
                axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

                if index_subplot <= 5

                    subplot(4,5,index_subplot+5)

                else

                   subplot(4,5,index_subplot+10) 

                end

                hold on
                box on
                grid on

                plot(sigma_xi_SR(variable_to_plot,:),'b','LineWidth',2)

                title(strcat('\sigma_{x',num2str(i),'} in proposed algorithm'))
                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';
                axis([0 Tsim 0 max(sigma_xi_SR(variable_to_plot,:))])

                index_subplot = index_subplot + 1;

            end

        end
    
    end

end

%%

if ~individual_plots
    
    if nx <= 10

        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        index_subplot = 1;

        if do_PD

            for i = 1:nx

                variable_to_plot = i;

                for t = 1:Tsim

                    upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_PD(variable_to_plot,t),sat_max_constraint);
                    lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_PD(variable_to_plot,t),sat_min_constraint);

                end

                if index_subplot <= 5

                    subplot(4,5,index_subplot)

                else

                   subplot(4,5,index_subplot+5) 

                end

                hold on
                box on
                grid on

                plot(x_PD(variable_to_plot,:),'b','LineWidth',2)
                plot(upper_bound,'b--','LineWidth',2)
                plot(lower_bound,'b--','LineWidth',2)

                title(strcat('x_{',num2str(i),'} in predifined route'))
                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';
                axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

                if index_subplot <= 5

                    subplot(4,5,index_subplot+5)

                else

                   subplot(4,5,index_subplot+10) 

                end

                hold on
                box on
                grid on

                plot(sigma_xi_PD(variable_to_plot,:),'b','LineWidth',2)

                title(strcat('\sigma_{x',num2str(i),'} in predifined route'))
                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';
                axis([0 Tsim 0 max(sigma_xi_PD(variable_to_plot,:))])

                index_subplot = index_subplot + 1;

            end

        end
        
    end

end

%%

if segments_used == 2

    if ~individual_plots

        if nx == 30

            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

            index_subplot = 1;

            if do_SR

                for i = 1:nx

                    variable_to_plot = i;

                    for t = 1:Tsim

                        upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_SR(variable_to_plot,t),sat_max_constraint);
                        lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_SR(variable_to_plot,t),sat_min_constraint);

                    end

                    subplot(5,6,index_subplot)

                    hold on
                    box on
                    grid on

                    plot(x_SR(variable_to_plot,:),'b','LineWidth',2)
                    plot(upper_bound,'b--','LineWidth',2)
                    plot(lower_bound,'b--','LineWidth',2)

                    title(strcat('x_{',num2str(i),'} in proposed algorithm'))
                    ax = gca;
                    ax.FontSize = 14;
                    ax.FontWeight = 'bold';
                    axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

                    index_subplot = index_subplot + 1;

                end

            end

        end

    end

end

%%

if segments_used == 2

    if ~individual_plots

        if nx == 30 

            figure
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

            index_subplot = 1;

            if do_PD

                for i = 1:nx

                    variable_to_plot = i;

                    for t = 1:Tsim

                        upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_PD(variable_to_plot,t),sat_max_constraint);
                        lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_PD(variable_to_plot,t),sat_min_constraint);

                    end

                    subplot(5,6,index_subplot)

                    hold on
                    box on
                    grid on

                    plot(x_PD(variable_to_plot,:),'b','LineWidth',2)
                    plot(upper_bound,'b--','LineWidth',2)
                    plot(lower_bound,'b--','LineWidth',2)

                    title(strcat('x_{',num2str(i),'} in predifined route'))
                    ax = gca;
                    ax.FontSize = 14;
                    ax.FontWeight = 'bold';
                    axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

                    index_subplot = index_subplot + 1;

                end

            end

        end

    end

end

%%

if segments_used == 1

    if ~individual_plots

        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        index_subplot = 1;

        if do_SR

            for i = 1:ns

                variable_to_plot = segments(i);

                for t = 1:Tsim

                    upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_SR(variable_to_plot,t),sat_max_constraint);
                    lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_SR(variable_to_plot,t),sat_min_constraint);

                end

                subplot(4,2,index_subplot)

                hold on
                box on
                grid on

                plot(x_real_SR(variable_to_plot,:),'b','LineWidth',2)
                plot(x_controller_SR(variable_to_plot,:),'g-.','LineWidth',2)
                plot(upper_bound,'b--','LineWidth',2)
                plot(lower_bound,'b--','LineWidth',2)

                %title(strcat('x_{',num2str(i),'} in proposed algorithm'))
                title(strcat('x_{',num2str(i),'}'))
                xlabel('time [step]')
                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';
                axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

                index_subplot = index_subplot + 1;

            end

        end
    
    end

end

%%

if segments_used == 1

    if ~individual_plots

        figure
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

        index_subplot = 1;

        if do_PD

            for i = 1:ns

                variable_to_plot = segments(i);

                for t = 1:Tsim

                    upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_PD(variable_to_plot,t),sat_max_constraint);
                    lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_PD(variable_to_plot,t),sat_min_constraint);
%                   upper_bound(t) = xmax(variable_to_plot);
%                   lower_bound(t) = xmin(variable_to_plot);

                end

                subplot(4,2,index_subplot)

                hold on
                box on
                grid on

                plot(x_real_PD(variable_to_plot,:),'b','LineWidth',2)
                plot(x_controller_PD(variable_to_plot,:),'g-.','LineWidth',2)
                plot(upper_bound,'b--','LineWidth',2)
                plot(lower_bound,'b--','LineWidth',2)

                title(strcat('x_{',num2str(i),'}'))
                xlabel('time [step]')
%                 ylabel(strcat('x_{',num2str(i),'} in predifined route')
                ax = gca;
                ax.FontSize = 14;
                ax.FontWeight = 'bold';
                axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

                index_subplot = index_subplot + 1;

            end

        end
        
    end

end


%% 

if individual_plots

    if do_SR

        for i = 1:nx
            
            figure

            variable_to_plot = i;

            for t = 1:Tsim

                upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_SR(variable_to_plot,t),sat_max_constraint);
                lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_SR(variable_to_plot,t),sat_min_constraint);

            end

            subplot(2,1,1)
            hold on
            box on
            grid on

            plot(x_real_SR(variable_to_plot,:),'b','LineWidth',2)
            plot(upper_bound,'b--','LineWidth',2)
            plot(lower_bound,'b--','LineWidth',2)

            title(strcat('x_{',num2str(i),'} in proposed algorithm'))
            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

            subplot(2,1,2)
            hold on
            box on
            grid on

            plot(sigma_xi_SR(variable_to_plot,:),'b','LineWidth',2)

            title(strcat('\sigma_{x',num2str(i),'} in proposed algorithm'))
            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([0 Tsim 0 max(sigma_xi_SR(variable_to_plot,:))])

        end
        
        for i = 1:nx
            
            figure

            variable_to_plot = i;

            for t = 1:Tsim

                upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_SR(variable_to_plot,t),sat_max_constraint);
                lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_SR(variable_to_plot,t),sat_min_constraint);

            end

            hold on
            box on
            grid on

            plot(x_real_SR(variable_to_plot,:),'b','LineWidth',2)
            plot(upper_bound,'b--','LineWidth',2)
            plot(lower_bound,'b--','LineWidth',2)

            title(strcat('x_{',num2str(i),'} in proposed algorithm'))
            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

        end

    end 

end

%%

if individual_plots

    if do_PD

        for i = 1:nx

            figure

            variable_to_plot = i;

            for t = 1:Tsim

                upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_PD(variable_to_plot,t),sat_max_constraint);
                lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_PD(variable_to_plot,t),sat_min_constraint);

            end

            subplot(2,1,1)
            hold on
            box on
            grid on

            plot(x_real_PD(variable_to_plot,:),'b','LineWidth',2)
            plot(upper_bound,'b--','LineWidth',2)
            plot(lower_bound,'b--','LineWidth',2)

            title(strcat('x_{',num2str(i),'} in predifined route'))
            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

            subplot(2,1,2)
            hold on
            box on
            grid on

            plot(sigma_xi_SR(variable_to_plot,:),'b','LineWidth',2)

            title(strcat('\sigma_{x',num2str(i),'} in predifined route'))
            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([0 Tsim 0 max(sigma_xi_SR(variable_to_plot,:))])

        end
        
        for i = 1:nx

            figure

            variable_to_plot = i;

            for t = 1:Tsim

                upper_bound(t) = max(xmax(variable_to_plot) - 3*sigma_xi_PD(variable_to_plot,t),sat_max_constraint);
                lower_bound(t) = min(xmin(variable_to_plot) + 3*sigma_xi_PD(variable_to_plot,t),sat_min_constraint);
 
            end
            
            hold on
            box on
            grid on

            plot(x_real_PD(variable_to_plot,:),'b','LineWidth',2)
            plot(upper_bound,'b--','LineWidth',2)
            plot(lower_bound,'b--','LineWidth',2)

            title(strcat('x_{',num2str(i),'} in predifined route'))
            ax = gca;
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            axis([0 Tsim xmin(variable_to_plot)-0.5 xmax(variable_to_plot)+0.5])

        end

    end    

end