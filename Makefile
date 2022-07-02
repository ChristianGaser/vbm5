OLDVERSION=v1.02
VERSION=v1.03

TARGET=/Users/gaser/spm/spm5/toolbox/vbm5

STARGET=141.35.200.101:/opt/lampp/htdocs/

FILES=cg_config_vbm.m cg_gwc_HMRF.m cg_preproc_write.m INSTALL.txt

ZIPFILE=vbm5_$(VERSION).zip

install: upgrade
	-@test ! -d ${TARGET} || rm -r ${TARGET}
	-@mkdir ${TARGET}
	-@cp ${FILES} ${TARGET}

upgrade:
	@for i in $(FILES); do sed -i "" -e "s/$(OLDVERSION)/$(VERSION)/g" $$i; done    

help:
	-@echo Available commands:
	-@echo install zip scp upgrade

zip: upgrade
	-@test ! -d vbm5 || rm -r vbm5
	-@cp -r ${TARGET} .
	-@zip ${ZIPFILE} -rm vbm5

scp: zip
	-@scp -pr ${ZIPFILE} ${STARGET}
