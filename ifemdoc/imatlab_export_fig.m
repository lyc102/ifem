function exported = imatlab_export_fig(exporter)
    % IMATLAB_EXPORT_FIG Set exporter or export figures for imatlab.
    %
    %   IMATLAB_EXPORT_FIG(exporter)
    %     where exporter is one of
    %       {'', 'fig2plotly', 'print-png', 'print-svg', 'print-jpeg'}
    %     sets the current exporter.
    %
    %   exported = IMATLAB_EXPORT_FIG
    %     orders the current figures by number, exports and closes them, and
    %     returns a cell array of exported filenames.

    persistent set_exporter
    if isempty(set_exporter)
        set_exporter = '';
    end
    valid_exporters = { ...
        '', 'fig2plotly', 'print-png', 'print-svg', 'print-jpeg'};

    if exist('exporter', 'var')
        if strcmp(exporter, '')
            set(0, 'defaultfigurevisible', 'on');
        else
            set(0, 'defaultfigurevisible', 'off');
        end
        if any(strcmp(exporter, valid_exporters))
            if strcmp(exporter, 'fig2plotly')
                version_delta = ...
                    str2double(strsplit(plotly_version, '.')) - [2, 2, 7];
                if version_delta(find(version_delta, 1)) < 0
                    error('imatlab:unsupportedPlotlyVersion', ...
                          'imatlab requires plotly>=2.2.7.')
                end
            end
            set_exporter = exporter;
        else
            error('imatlab:invalidExporter', ...
                  ['known exporters are ', ...
                   strjoin(cellfun(@(c) ['''', c, ''''], valid_exporters, ...
                           'UniformOutput', false), ', ')]);
        end
    elseif strcmp(set_exporter, '')
        exported = {};
    else
        children = get(0, 'children');
        [~, idx] = sort([children.Number]);
        children = children(idx);
        exported = cell(1, numel(children));
        for i = 1:numel(children)
            child = children(i);
            name = tempname('.');
            if strcmp(set_exporter, 'fig2plotly')
                exported{i} = [name, '.html'];
                try
                    fig2plotly(child, 'filename', name, ...
                               'offline', true, 'open', false);
                catch me
                    warning('fig2plotly failed to export a figure');
                    rethrow(me);
                end
            else
                ihc = child.InvertHardcopy;
                child.InvertHardcopy = 'off';  % Respect user background.
                ext = set_exporter(1+numel('print-'):end);
                exported{i} = [name, '.', ext];
                % Use screen resolution.
                print(child, exported{i}, ['-d', ext], '-r0');
                child.InvertHardcopy = ihc;
            end
            close(child);
        end
    end
end