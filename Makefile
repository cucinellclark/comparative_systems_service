TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime

APP_SERVICE = app_service

WRAP_PYTHON_TOOL = wrap_python3
WRAP_PYTHON_SCRIPT = bash $(TOOLS_DIR)/$(WRAP_PYTHON3_TOOL).sh

SRC_PYTHON = $(wildcard scripts/*.py)

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

all: bin

bin: $(BIN_PYTHON) $(BIN_SERVICE_PERL)

deploy: deploy-client
deploy-all: deploy-client
deploy-client: deploy-scripts deploy-libs

deploy-service: deploy-libs deploy-scripts deploy-service-scripts deploy-specs

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs:

clean:

include $(TOP_DIR)/tools/Makefile.common.rules
