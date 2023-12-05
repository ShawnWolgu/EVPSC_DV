local dap = require('dap')
dap.adapters.codelldb = {
    type = 'server';
    port = '${port}';
    executable = {
        command = vim.fn.exepath('codelldb'),
		args = { '--port', '${port}'},
    };
}
dap.configurations.cpp = {
    {
        type = 'codelldb';
        request = "launch";
        name = 'Debug';
        program = '${workspaceFolder}/build/EVPSC';
        args = {};
        cwd = '${workspaceFolder}/debug/';
        terminal = 'integrated';
        console = 'integratedTerminal';
    }
}

