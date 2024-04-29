
# export ALIBUILD_WORK_DIR, FRAGMENTATION_DIR wherever you call this from
# eval "`alienv shell-helper`" wherever you call this from
# equip with conda loader if needed

# Change prompt to show environment name
if [ x$O2_ROOT != x ]
then
export PS1='%T %F{208}[O2 ninja] %F{magenta}%~%F{105}${vcs_info_msg_0_}%F{208} %#> %f'
else
export PS1='%T %F{magenta}%~%F{105}${vcs_info_msg_0_}%f %# '
fi
export CLICOLOR=1
export LSCOLORS=gxgxBxDxCxEgEdxbxgxcxd

# Version control settings
autoload -Uz compinit && compinit
autoload -Uz vcs_info
precmd() { vcs_info }
zstyle ':vcs_info:git*' formats ' (%b)'
setopt PROMPT_SUBST

source .aliases
source .zsh_aliases
