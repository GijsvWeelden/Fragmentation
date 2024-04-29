" Set line numbers
set number
" Make line numbers not selectable by mouse
" set mouse=n
syntax on
" Indent settings
filetype plugin indent on
set smartindent
set expandtab
set tabstop=2
set shiftwidth=2
" autocmd BufRead,BufWritePre *.sh normal gg=G
" autocmd BufRead,BufWritePre *.C normal gg=G
" autocmd BufRead,BufWritePre *.cxx normal gg=G
" autocmd BufRead,BufWritePre *.h normal gg=G
" autocmd BufRead,BufWritePre *.hh normal gg=G
" autocmd BufRead,BufWritePre *.cc normal gg=G
" Automatically strip trainling whitespace when a file is saved
:autocmd BufWritePre * :%s/\s\+$//e
" Show when going beyond 100 characters on a line
highlight ColorColumn ctermbg=LightGray
call matchadd('ColorColumn', '\%101v', 100)
highlight CursorLine cterm=NONE ctermbg=NONE ctermfg=NONE guibg=NONE guifg=NONE
set cursorline
" Highlight search hits
set hlsearch
" Turn off search highlight when entering insert mode
:nnoremap i :noh<cr>i
