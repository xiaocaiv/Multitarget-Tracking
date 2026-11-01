function params = apply_scenario(params, name)
%APPLY_SCENARIO 把预置场景写入 params (目标初始状态 / 机动段 / 变速段)
%   params = APPLY_SCENARIO(params, name)
%
%   name 取自 scenario_defs.m 的 .name; 写入:
%     params.target1_init, params.target2_init
%     params.target_turns, params.target_accel
%     params.maneuver_enable (场景含机动/变速段时为 true)
%     params.scenario, params.scenario_label

S = scenario_defs();
idx = find(strcmpi({S.name}, name), 1);
if isempty(idx)
    error('apply_scenario:unknown', '未知场景 "%s"; 可选: %s', ...
        name, strjoin({S.name}, ', '));
end
sc = S(idx);

params.scenario = sc.name;
params.scenario_label = sc.label;
params.target1_init = sc.targets{1};
params.target2_init = sc.targets{2};
params.target_turns = sc.turns;
params.target_accel = sc.accel;

has_turn = any(~cellfun(@isempty, sc.turns));
has_accel = any(~cellfun(@isempty, sc.accel));
params.maneuver_enable = has_turn || has_accel;
end
