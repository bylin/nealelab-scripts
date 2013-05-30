" Automatically cd into the directory the file is in
autocmd BufEnter * execute "chdir ".escape(expand("%:p:h"), ' ')
set showcmd
set foldmethod=marker
set softtabstop=3
set autoindent
set shiftwidth=3
set mouse=a
set nu
set incsearch
set hlsearch
set modeline
set ls=2
